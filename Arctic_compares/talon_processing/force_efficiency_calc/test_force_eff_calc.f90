program test_force_eff_calc
!
! NAME:
!   test_force_eff_calc
!
! PURPOSE:
! 
! CALLS:
!   Modules:
!     - hdf5
!   Subroutines:
!     - check_bad_rows
!
! MODIFICATIONS:
!   Blake Sorenson <blake.sorenson@und.edu>     - 2024/02/29:
!     Written
!
!  ############################################################################

  use hdf5
  use comp_grid_vars, only: clear_arrays, i_bad_list, &
    find_distance_between_points, allocate_out_arrays, &
    pixel_in_box, &
    MODIS_CH7_data, MODIS_CLD_data, MODIS_COD_data, &
    CERES_SWF_data, OMI_AI_data, &
    OMI_SZA_data, NSIDC_data, OMI_AI_dims
  use force_eff_vars, only: clear_force_eff_arrays, allocate_force_out_arrays, &
    sza_bins, ice_bins, cod_bins, sza_edges, ice_edges, cod_edges, &
    sumx_values, sumy_values, sumx2_values, sumy2_values, sumxy_values, &
    count_values, regress_slopes, regress_intercepts, regress_slope_error, &
    regress_r2, len_sza, len_ice, len_cld, len_cod

  implicit none

  real                  :: min_ai
  real                  :: max_ai
  real                  :: min_sza
  real                  :: max_sza
  real                  :: delta_sza
  real                  :: min_ice
  real                  :: max_ice
  real                  :: delta_ice
  real                  :: min_cod
  real                  :: max_cod
  real                  :: delta_cod

  integer               :: mm
  integer               :: nn
  integer               :: ii
  integer               :: jj

  integer               :: sza_idx
  integer               :: ice_idx
  integer               :: cld_idx
  integer               :: cod_idx

  integer               :: check_sza_idx
  integer               :: check_ice_idx
  integer               :: check_cld_idx
  integer               :: check_cod_idx

  integer               :: rank
  integer               :: i_file_length
  integer(hsize_t)      :: file_id
  integer(hsize_t)      :: out_file_id
  integer(hsize_t)      :: dspace_id
  integer(hsize_t)      :: dset_id
  integer(hsize_t), dimension(:)                   :: regress_dims(4)

  !!#!real, dimension(:), allocatable                  :: sza_bins
  !!#!real, dimension(:), allocatable                  :: ice_bins
  !!#!real, dimension(:), allocatable                  :: cod_bins

  !!#!real, dimension(:), allocatable                  :: sza_edges
  !!#!real, dimension(:), allocatable                  :: ice_edges
  !!#!real, dimension(:), allocatable                  :: cod_edges

  !!#!real, dimension(:,:,:,:), allocatable                  :: sumx_values
  !!#!real, dimension(:,:,:,:), allocatable                  :: sumy_values
  !!#!real, dimension(:,:,:,:), allocatable                  :: sumx2_values
  !!#!real, dimension(:,:,:,:), allocatable                  :: sumy2_values
  !!#!real, dimension(:,:,:,:), allocatable                  :: sumxy_values
  !!#!integer, dimension(:,:,:,:), allocatable               :: count_values

  real                                                   :: meanx_values
  real                                                   :: meany_values

  !!#!real, dimension(:,:,:), allocatable                  :: sum_pred_err_square
  !!#!real, dimension(:,:,:), allocatable                  :: sum_x_pert_square

  real                                                   :: s_xy
  real                                                   :: s_xx
  real                                                   :: sse 
  real                                                   :: sst 

  !!#!real, dimension(:,:,:,:), allocatable                  :: regress_slopes
  !!#!real, dimension(:,:,:,:), allocatable                  :: regress_intercepts
  !!#!real, dimension(:,:,:,:), allocatable                  :: regress_slope_error
  !!#!real, dimension(:,:,:,:), allocatable                  :: regress_r2

  real, dimension(4)    :: ch7_thresh
  real                  :: local_thresh

  real                  :: avg_ai

  integer               :: num_misclass

  integer               :: min_ice_idx
  integer               :: min_sza_idx
  integer               :: min_cld_idx
  integer               :: min_cod_idx
  integer               :: max_ice_idx
  integer               :: max_sza_idx
  integer               :: max_cld_idx
  integer               :: max_cod_idx

  real                  :: y_hat

  integer               :: minii
  integer               :: maxii
  integer               :: minjj
  integer               :: maxjj

  integer               :: arg_count
  integer               :: io7
  integer               :: io8
  integer               :: istatus
  integer               :: error
  character(len = 255)  :: file_name_file
  character(len = 12)   :: single_file_date
  character(len = 255)  :: total_file_name
  character(len = 255)  :: out_file_name
  
  ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  arg_count = command_argument_count()
  if(arg_count /= 1) then
    write(*,*) 'SYNTAX: ./force_eff_calc_exec file_name_file'
    write(*,*) '         file_name_file consists of lines of YYYYMMDDHHMM'
    return
  endif

  call get_command_argument(1, file_name_file)

  ! Initialize the HDF5 interface
  ! -----------------------------
  call h5open_f(error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: Could not open HDF5 library'
    return
  endif
  write(*,*) 'Interface opened'

  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !
  ! Step 1:  Determine how many bins for each variable are needed
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  min_ice_idx = 1000
  min_sza_idx = 1000
  min_cld_idx = 1000
  min_cod_idx = 1000
  max_ice_idx = -1000
  max_sza_idx = -1000
  max_cld_idx = -1000
  max_cod_idx = -1000

  ! Specify the bin ranges and sizes
  ! --------------------------------
  min_ai   = -2.0
  max_ai   = 12.0

  min_sza   = 40.0
  max_sza   = 85.0
  delta_sza = 5.0

  min_ice   = 0.0
  max_ice   = 105.0
  delta_ice = 5.0

  min_cod   = 0.01
  max_cod   = 40.01
  delta_cod = 2.0

  ! Define bin centers and bin edges for each of these
  ! --------------------------------------------------
  len_sza = (max_sza - min_sza + delta_sza) / delta_sza
  !len_ice = (max_ice - min_ice + delta_ice) / delta_ice
  len_cod = 8
  !len_cod = ((max_cod - min_cod + delta_cod) / delta_cod) + 1
  len_ice = 4
  len_cld = 2

  allocate(sza_bins(len_sza))
  allocate(ice_bins(len_ice))
  allocate(cod_bins(len_cod))
  !allocate(cld_bins(len_ice))
  !allocate(cld_bins(len_cld))

  allocate(sza_edges(len_sza + 1))
  allocate(ice_edges(len_ice + 1))
  allocate(cod_edges(len_cod + 1))
  !allocate(ch7_edges(len_ch7 + 1))

  do ii = 1, len_sza
    sza_bins(ii) = min_sza + delta_sza*(ii-1)
    write(*,*) 'SZA bin',ii,sza_bins(ii)
  enddo

  ! Test using the Feng and Christopher COD bins
  ! --------------------------------------------
  cod_bins(1) = 0. 
  cod_bins(2) = 1.
  cod_bins(3) = 3.
  cod_bins(4) = 5.
  cod_bins(5) = 7.
  cod_bins(6) = 10.
  cod_bins(7) = 16.
  cod_bins(8) = 26.

  cod_edges(1) = 0.
  cod_edges(2) = 0.00001
  cod_edges(3) = 2.0
  cod_edges(4) = 4.0
  cod_edges(5) = 6.0
  cod_edges(6) = 8.0
  cod_edges(7) = 12.0
  cod_edges(8) = 20.
  cod_edges(9) = 70.

  !write(*,*) 'COD bin',1,cod_bins(1)
  !do ii = 2, len_cod
  !  !write(*,*) ii, len_cod, min_cod + delta_cod*(ii-2)
  !  cod_bins(ii) = min_cod + delta_cod*(ii-2)
  !  write(*,*) 'COD bin',ii,cod_bins(ii)
  !enddo

  write(*,*) "HERE1",cod_bins(1:)
  write(*,*) "HERE2",cod_bins(2:)
 
  ice_bins = [10, 90, 105, 50]

  do ii = 1, len_ice
    !!#!ice_bins(ii) = min_ice + delta_ice*(ii-1)
    write(*,*) 'ICE bin',ii,ice_bins(ii)
  enddo


  !do ii = 1, len_cld
  !  cld_bins(ii) = min_cld + delta_cld*(ii-1)
  !enddo

  ice_edges(1) = 0.
  ice_edges(2) = 20.
  ice_edges(3) = 80.
  ice_edges(4) = 100.
  ice_edges(5) = 254.
  do ii = 1, len_sza + 1 
    sza_edges(ii) = min_sza + delta_sza*((ii-1) - 0.5)
  enddo
 
  !!#!do ii = 1, len_ice + 1 
  !!#!  ice_edges(ii) = min_ice + delta_ice*((ii-1) - 0.5)
  !!#!enddo

  ! Allocate arrays to hold the 2.1 μm thresholds for splitting
  ! between clear and cloudy conditions
  ! The indices follow the ice type:
  !   1: ocean
  !   2: ice
  !   3: land
  !   4: mix
  ! -----------------------------------------------------------
  ch7_thresh = [0.03,0.035,0.07, 0.03] 

  ! Use these bins to allocate and initialize arrays to hold
  ! all of the following parameters:
  allocate(sumx_values(len_cod, len_cld, len_sza, len_ice))
  allocate(sumy_values(len_cod, len_cld, len_sza, len_ice))
  allocate(sumx2_values(len_cod, len_cld, len_sza, len_ice))
  allocate(sumy2_values(len_cod, len_cld, len_sza, len_ice))
  allocate(sumxy_values(len_cod, len_cld, len_sza, len_ice))
  allocate(count_values(len_cod, len_cld, len_sza, len_ice))

  sumx_values(:,:,:,:) = 0.
  sumy_values(:,:,:,:) = 0.
  sumx2_values(:,:,:,:) = 0.
  sumy2_values(:,:,:,:) = 0.
  sumxy_values(:,:,:,:) = 0.
  count_values(:,:,:,:) = 0


  ! This can be calculated after the first loop and  before
  ! the second. This is needed for determining the slope
  ! standard error
  !allocate(meanx_values(len_cod, len_cld, len_sza, len_ice))
  !allocate(meany_values(len_cod, len_cld, len_sza, len_ice))
  !allocate(s_xy(len_cod, len_cld, len_sza, len_ice))
  !allocate(s_xx(len_cod, len_cld, len_sza, len_ice))
  !allocate(sse(len_cod, len_cld, len_sza, len_ice))
  !allocate(sst(len_cod, len_cld, len_sza, len_ice))

  !!#!! These are used in the determining of the slope standard error
  !!#!! sum_pred_err_square: Σ(yi - y_hati)**2.
  !!#!allocate(sum_pred_err_square(len_cld, len_sza, len_ice))
  !!#!! sum_x_pert: Σ(xi - x_mean)**2.
  !!#!allocate(sum_x_pert_square(len_cld, len_sza, len_ice))

  !meanx_values(:,:,:,:) = 0.
  !meany_values(:,:,:,:) = 0.
  !s_xy(:,:,:,:)         = 0.
  !s_xx(:,:,:,:)         = 0.
  !sse(:,:,:,:)          = 0.
  !sst(:,:,:,:)          = 0.
  !!#!sum_pred_err_square(:,:,:) = 0.
  !!#!sum_x_pert_square(:,:,:) = 0.

  allocate(regress_slopes(len_cod, len_cld, len_sza, len_ice))
  allocate(regress_intercepts(len_cod, len_cld, len_sza, len_ice))
  allocate(regress_slope_error(len_cod, len_cld, len_sza, len_ice))
  allocate(regress_r2(len_cod, len_cld, len_sza, len_ice))

  regress_slopes(:,:,:,:) = 0.
  regress_intercepts(:,:,:,:) = 0.
  regress_slope_error(:,:,:,:) = 0.
  regress_r2(:,:,:,:) = 0.

  !if(  (nn == check_cld_idx) .and. (jj == check_sza_idx) .and. &
  !     (ii == check_ice_idx) .and. (mm == check_cod_idx) ) then

  check_cod_idx = 1
  check_cld_idx = 1
  check_sza_idx = 5
  check_ice_idx = 4
  !check_sza_idx = 7 !- SZA = 65
  !check_ice_idx = 1 !- ICE = 10 (ocean)
  !check_cld_idx = 1 !- CLD = clear
  !check_sza_idx = 4 !- SZA = 50
  !check_ice_idx = 1 !- ICE = 10 (ocean)
  !check_cld_idx = 2 !- CLD = cloudy

  ! = = = = = = = = = = = = = = = = = = = = = = = = = =
  !
  ! Step 2: Loop over all the colocated data products and calculate, for each
  !         surface type and SZA, the following
  !         (y = ceres_swf, x = omi_ai)
  !   
  !
  !     Use the MODIS CLD cloudmask, surrounding AI, surface type,
  !       and others to determine if this pixel is actually cloudy or
  !       if it is misclassified 
  !
  !     - Σx      : sum of all the OMI AI in these conditions
  !     - Σy      : sum of all the CERES SWF in these conditions
  !     - Σx2     : sum of the squares of all the OMI AI in these conditions
  !     - Σxy     : sum of each CERES SWF multiplied by each OMI AI
  !     - (Σx)^2  : square of the sum of all the OMI AI in these conditions
  !     - nn      : number of pixels in this classification type
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = = = =

  ! Open the file name file
  io7 = 2007
  open(io7, file = 'test_out_file.txt', iostat = istatus)
  if(istatus > 0) then
    write(*,*) 'ERROR: Problem opening test_out_file.txt'
  endif
  write(io7,'(a9)') 'OMI,CERES' 

  ! Open the file name file
  io8 = 2023
  open(io8, file = trim(file_name_file), iostat = istatus)
  if(istatus > 0) then
    write(*,*) 'ERROR: Problem reading '//trim(file_name_file)
    return
  else

    ! Find the length of the input file
    ! ---------------------------------
    i_file_length = 0
    do
      read(io8, *, iostat = istatus)
      if(istatus < 0) then
        exit
      else 
        i_file_length = i_file_length + 1
      endif
    enddo
   
    write(*,*) 'Length of the input file:', i_file_length
 
    rewind(io8) 


    num_misclass = 0
    file_loop: do nn = 1, i_file_length

      ! Open the current colocated file
      read(io8, '(A)', iostat = istatus) single_file_date
      if(istatus > 0) then
        write(*,*) 'ERROR: Problem reading line from '//trim(file_name_file)
        cycle file_loop
      else
        !!#!write(*,*) nn, single_file_date

        total_file_name = '/home/blake.sorenson/OMI/arctic_comp/comp_data/'&
            //'colocated_subset_'//single_file_date//'.hdf5'

        ! Read data from this file
        ! ------------------------
        call h5fopen_f(total_file_name, H5F_ACC_RDWR_F, file_id, error)
        if(error /= 0) then
          write(*,*) 'FATAL ERROR: Could not open file '//total_file_name
          cycle file_loop
        endif

        call read_coloc_OMI_AI(file_id) 
        call read_coloc_OMI_SZA(file_id) 
        call read_coloc_NSIDC_ICE(file_id) 
        call read_coloc_MODIS_COD(file_id) 
        call read_coloc_CERES_SWF(file_id) 
        call read_coloc_MODIS_CLD(file_id) 
        call read_coloc_MODIS_CH7(file_id) 

        !!#!minii = 9
        !!#!maxii = 11
        !!#!minjj = 1
        !!#!maxjj = 2
        !!#!write(*,*) minii, maxii, minjj, maxjj

        !!#!write(*,*) OMI_AI_data(minjj:maxjj, minii:maxii)
        !!#!write(*,*) sum(OMI_AI_data(minjj:maxjj, minii:maxii)), &
        !!#!    size(OMI_AI_data(minjj:maxjj, minii:maxii)), &
        !!#!    sum(OMI_AI_data(minjj:maxjj, minii:maxii)) /  &
        !!#!    size(OMI_AI_data(minjj:maxjj, minii:maxii))

        omi_loop1: do ii = 1, OMI_AI_dims(2)

          ! Check if the current pixel is along an edge
          if(ii == 1) then
            minii = 1
            maxii = ii + 1
          else if(ii == OMI_AI_dims(2)) then
            minii = ii - 1
            maxii = ii
          else
            minii = ii - 1
            maxii = ii + 1
          endif

          omi_loop2: do jj = 1, OMI_AI_dims(1)

            ! Make sure missing AI values are not included
            ! --------------------------------------------
            if( ( (OMI_AI_data(jj,ii) > 2) .and. &
                  (OMI_AI_data(jj,ii) < 12) ) .and. &
                ( (CERES_SWF_data(jj,ii) > -200.) .and. &
                  (CERES_SWF_data(jj,ii) < 3000) ) .and. &
                ( MODIS_COD_data(jj,ii) /= -999.) .and. &
                ( NSIDC_data(jj,ii) /= 253.) ) then

              ! Identify the COD value
              ! ----------------------
              if(MODIS_COD_data(jj,ii) == 0.) then
                cod_idx = 1
              else

                !!#!! Add the '2:' and '+ 1 ' to this index check to ensure that
                !!#!! none of this data goes into index 0, which is reserved
                !!#!! only for COD exactly equal to 0.
                !!#!!cod_idx = minloc(abs(cod_bins - MODIS_COD_data(jj,ii)), dim = 1)
                !!#!cod_idx = minloc(abs(cod_bins(2:) - MODIS_COD_data(jj,ii)), dim = 1)  + 1


                cod_idx = -1
                do mm = 2, len_cod 
                  if(  (MODIS_COD_data(jj,ii) >= cod_edges(mm)) .and. &
                       (MODIS_COD_data(jj,ii) < cod_edges(mm+1)) ) then
                    cod_idx = mm
                  endif
                enddo
       
                ! Cycle the OMI loop to the next observation if the current
                ! COD data does not fall into any of the bins
                if( cod_idx == -1) then
                  cycle omi_loop2
                endif

              endif

              ! Identify the SZA and ICE values
              ! -------------------------------
              sza_idx = minloc(abs(sza_bins - OMI_SZA_data(jj,ii)), dim = 1) 
              !ice_idx = minloc(abs(ice_bins - NSIDC_data(jj,ii)), dim = 1) 

              if( (NSIDC_data(jj,ii) >= 0.) .and. (NSIDC_data(jj,ii) <= 20.) ) then
                ice_idx = 1
              else if( (NSIDC_data(jj,ii) > 20.) .and. (NSIDC_data(jj,ii) < 80.)) then
                ice_idx = 4
              else if( (NSIDC_data(jj,ii) > 80.) .and. (NSIDC_data(jj,ii) <= 100.)) then
                ice_idx = 2
              else if( NSIDC_data(jj,ii) == 254.) then
                ice_idx = 3
              else
                !!#!write(*,*) 'INVALID ICE VALUE. Skipping', NSIDC_data(jj,ii)
                cycle omi_loop2
              endif
                

              ! Apply the cloud mask check to determine if this pixel is cloudy or 
              ! clear
              !
              ! 3 == clear
              ! 2 == probably clear
              ! 1 == probably cloudy
              ! 0 == cloudy

              ! Am using this to figure out when the MODIS cloud mask classifies bins 
              ! as "cloudy" when the COD is exactly 0.
              !if( (MODIS_CLD_data(jj,ii) == 0.) .and. (MODIS_COD_data(jj,ii) == 0.)) then
              !  write(*,*) single_file_date
              !endif

                   
              ! If the bin is classified as clear, continue 
              if(MODIS_CLD_data(jj,ii) == 3.) then
                ! set cloud_index to 0 and continue
                cld_idx = 1    
                
              else
                ! Calculate the average of the AI values around this pixel
                ! --------------------------------------------------------
                ! Check if the current pixel is along an edge
                if(jj == 1) then
                  minjj = 1
                  maxjj = jj + 1
                else if(jj == OMI_AI_dims(1)) then
                  minjj = jj - 1
                  maxjj = jj
                else
                  minjj = jj - 1
                  maxjj = jj + 1
                endif


                avg_ai = sum(OMI_AI_data(minjj:maxjj, minii:maxii)) /  &
                        size(OMI_AI_data(minjj:maxjj, minii:maxii))
              
                ! Check the AI threshold to see if this is a significantly smoky 
                ! region
                if( (OMI_AI_data(jj,ii) > 1.5 ) .and. ( avg_ai > 1.5 )) then

                  ! Use the NSIDC_ICE value to figure out which threshold
                  ! value to use
                  ! -----------------------------------------------------
                  if(NSIDC_data(jj,ii) <= 20.) then
                    local_thresh = ch7_thresh(1)
                  else if((NSIDC_data(jj,ii) > 20.) .and. &
                          (NSIDC_data(jj,ii) < 80)) then
                    local_thresh = ch7_thresh(4)
                  else if((NSIDC_data(jj,ii) >= 80.) .and. &
                          (NSIDC_data(jj,ii) <= 100)) then
                    local_thresh = ch7_thresh(2)
                  else if(NSIDC_data(jj,ii) == 254.) then
                    local_thresh = ch7_thresh(3)
                  endif
                   

                  ! Check the MODIS 2.1 μm reflectance to see if this is a 
                  ! misclassification
                  if((MODIS_CH7_data(jj,ii) < local_thresh) .and. &
                     (MODIS_CH7_data(jj,ii) /= -999.)) then
                    ! MISCLASSIFICATION of cloudy as clear. Set cloud_index to 
                    ! 0 and continue
                    cld_idx = 1

                    ! Also, set the COD index to 1 as well, since the COD 
                    ! could also be contaminated by the misclassification
                    ! ---------------------------------------------------
                    cod_idx = 1
                 
                    num_misclass = num_misclass + 1
 
                  else
                    ! Correct class. of cloudy as cloudy. Set cloud_index to 1 
                    ! and continue
                    cld_idx = 2
               
                  endif    
              
                else 
                  ! The smoke here is not very high, so trust the initial cloud mask
                  ! and set cloud index to 'cloudy' 1
                  cld_idx = 2    
              
                endif
              endif

              ! Add a check to see if the cloud classification is "cloudy" and
              ! the COD is exactly 0. If it is, set the cloud flag to "clear"
              ! --------------------------------------------------------------
              if( (cld_idx == 2) .and. (MODIS_COD_data(jj,ii) == 0.)) then
                cld_idx = 1
              endif     

              if(ice_idx < min_ice_idx) then
                min_ice_idx = ice_idx
              else if(ice_idx > max_ice_idx) then
                max_ice_idx = ice_idx
              endif

              if(sza_idx < min_sza_idx) then
                min_sza_idx = sza_idx
              else if(sza_idx > max_sza_idx) then
                max_sza_idx = sza_idx
              endif

              if(cld_idx < min_cld_idx) then
                min_cld_idx = cld_idx
              else if(cld_idx > max_cld_idx) then
                max_cld_idx = cld_idx
              endif

              if(cod_idx < min_cod_idx) then
                min_cod_idx = cod_idx
              else if(cod_idx > max_cod_idx) then
                max_cod_idx = cod_idx
              endif


              !!#!if((cod_idx == 2) .and. (MODIS_CLD_data(jj,ii) == 3.)) then
              !!#!  write(*,*) 'COD index',cod_idx,'COD value',MODIS_COD_data(jj,ii)
              !!#!endif
              if(  (cld_idx == check_cld_idx) .and. (sza_idx == check_sza_idx) .and. &
                   (ice_idx == check_ice_idx) .and. (cod_idx == check_cod_idx) ) then
                !write(*,*) OMI_AI_data(jj,ii),',',CERES_SWF_data(jj,ii)
                write(*,*) 'OMI AI = ', OMI_AI_data(jj,ii), ' CERES SWF = ', &
                    CERES_SWF_data(jj,ii)
                write(io7,'(f5.3,1a,f6.2)') OMI_AI_data(jj,ii),',',CERES_SWF_data(jj,ii)
              endif

              sumx_values(cod_idx,cld_idx, sza_idx, ice_idx) = &
                sumx_values(cod_idx,cld_idx, sza_idx, ice_idx) + OMI_AI_data(jj,ii)
              sumy_values(cod_idx,cld_idx, sza_idx, ice_idx) = &
                sumy_values(cod_idx,cld_idx, sza_idx, ice_idx) + CERES_SWF_data(jj,ii)
              sumx2_values(cod_idx,cld_idx, sza_idx, ice_idx) = &
                sumx2_values(cod_idx,cld_idx, sza_idx, ice_idx) + OMI_AI_data(jj,ii)**2.
              sumy2_values(cod_idx,cld_idx, sza_idx, ice_idx) = &
                sumy2_values(cod_idx,cld_idx, sza_idx, ice_idx) + CERES_SWF_data(jj,ii)**2.
              sumxy_values(cod_idx,cld_idx, sza_idx, ice_idx) = &
                sumxy_values(cod_idx,cld_idx, sza_idx, ice_idx) + &
                (OMI_AI_data(jj,ii) * CERES_SWF_data(jj,ii))
              count_values(cod_idx,cld_idx, sza_idx, ice_idx) = &
                count_values(cod_idx,cld_idx, sza_idx, ice_idx) + 1

            endif
          enddo omi_loop2
        enddo omi_loop1

        ! Deallocate all arrays for the next pass
        ! ---------------------------------------
        call clear_arrays

        ! Close file
        ! ----------
        call h5fclose_f(file_id, error)
        if(error /= 0) then
          write(*,*) 'FATAL ERROR: Could not close file'
          return
        endif

      endif
    enddo file_loop
  endif

  close(io7)

  write(*,*) 'Shape of sumx_values', shape(sumx_values)
  write(*,*) 'Min/max ice idx', min_ice_idx, max_ice_idx
  write(*,*) 'Min/max sza idx', min_sza_idx, max_sza_idx
  write(*,*) 'Min/max cld idx', min_cld_idx, max_cld_idx
  write(*,*) 'Min/max cod idx', min_cod_idx, max_cod_idx

  write(*,*) 'Num misclass:', num_misclass

  ! = = = = = = = = = = = = = = = = = = = = = = = = = =
  !
  ! Step 3: For each of these conditions, use the accumulated parameters to
  !         calculate the slope and intercept of the regression lines
  !         
  ! = = = = = = = = = = = = = = = = = = = = = = = = = =
  do ii = 1, len_ice
    do jj = 1, len_sza
      do nn = 1, len_cld
        do mm = 1, len_cod

          !!#!if(count_values(nn,jj,ii) > 1) then
          !!#!  regress_intercepts(nn,jj,ii) = ((sumy_values(nn,jj,ii) * &
          !!#!      (sumx_values(nn,jj,ii) **2.)) - (sumx_values(nn,jj,ii) * &
          !!#!       sumxy_values(nn,jj,ii))) / &
          !!#!      (count_values(nn,jj,ii) * sumx2_values(nn,jj,ii) - &
          !!#!       (sumx_values(nn,jj,ii)**2.))

          !!#!  regress_slopes(nn,jj,ii) = ((count_values(nn,jj,ii) * &
          !!#!      sumxy_values(nn,jj,ii)) - (sumx_values(nn,jj,ii) - &
          !!#!      sumy_values(nn,jj,ii))) / &
          !!#!      (count_values(nn,jj,ii) * sumx2_values(nn,jj,ii) - &
          !!#!       (sumx_values(nn,jj,ii)**2.))

          !!#!  write(*,*) nn,jj,ii, count_values(nn,jj,ii), &
          !!#!    regress_slopes(nn,jj,ii), regress_intercepts(nn,jj,ii)
          !!#!endif

          ! Calculate the mean x (AI) value here per bin
          ! --------------------------------------------
          if(count_values(mm,nn,jj,ii) > 1) then
            meanx_values = sumx_values(mm,nn,jj,ii) / count_values(mm,nn,jj,ii)
            meany_values = sumy_values(mm,nn,jj,ii) / count_values(mm,nn,jj,ii)


            s_xy = sumxy_values(mm,nn,jj,ii) - &
                ( (sumx_values(mm,nn,jj,ii) * sumy_values(mm,nn,jj,ii)) / &
                count_values(mm,nn,jj,ii))
            s_xx = sumx2_values(mm,nn,jj,ii) - &
                ( (sumx_values(mm,nn,jj,ii)**2.) / count_values(mm,nn,jj,ii) ) 

            regress_slopes(mm,nn,jj,ii) = s_xy / s_xx

            regress_intercepts(mm,nn,jj,ii) = meany_values - &
                regress_slopes(mm,nn,jj,ii) * meanx_values

            sse = sumy2_values(mm,nn,jj,ii) - &
                regress_intercepts(mm,nn,jj,ii) * sumy_values(mm,nn,jj,ii) - &
                regress_slopes(mm,nn,jj,ii) * sumxy_values(mm,nn,jj,ii)

            sst = sumy2_values(mm,nn,jj,ii) - &
                (sumy_values(mm,nn,jj,ii)**2.) / count_values(mm,nn,jj,ii)

            regress_r2(mm,nn,jj,ii) = 1. - (sse / sst)
 
            regress_slope_error(mm,nn,jj,ii) = ( (1. / &
                (count_values(mm,nn,jj,ii) - 2)) * &
                ( sse / s_xx  )  ) ** 0.5

            if(  (nn == check_cld_idx) .and. (jj == check_sza_idx) .and. &
                 (ii == check_ice_idx) .and. (mm == check_cod_idx) ) then
              write(*,*) 'CHECKING THE FOLLOWING'
            endif
            write(*,'(4(i4),f7.2,1x,i3,1x,f5.1,1x,f6.2,i7,1x,f9.2,2x,f7.2,2x,f7.2)') &
                mm, nn, jj, ii, &
                cod_bins(mm), nn - 1, sza_bins(jj), ice_bins(ii),&
                count_values(mm,nn,jj,ii), regress_slopes(mm,nn,jj,ii), &
                regress_slope_error(mm,nn,jj,ii), regress_r2(mm,nn,jj,ii)

          endif
        enddo
      enddo
    enddo
  enddo

  !!#!regress_intercept = ((Σy * (Σx)**2. ) - (Σx * Σxy)) / &
  !!#!                        nn * Σx2 - (Σx)**2.

  !!#!regress_slope     = (nn * Σxy - (Σx * Σy)) / &
  !!#!                     nn * Σx2 - (Σx)**2.

  !!#!rewind(io8)
 
  !!#!!!#!! = = = = = = = = = = = = = = = = = = = = = = = = = =
  !!#!!!#!!
  !!#!!!#!! Step 4: Loop back over all the colocated data products \
  !!#!!!#!!         and calculate, for each
  !!#!!!#!!         surface type and SZA, the following
  !!#!!!#!!         (y = ceres_swf, x = omi_ai)
  !!#!!!#!!
  !!#!!!#!!
  !!#!!!#!!         Need the following:
  !!#!!!#!!
  !!#!!!#!!       y_hat              : predicted SWF value given the regression and AI
  !!#!!!#!!       Σ(yi - y_hati)**2.
  !!#!!!#!!       Σ(xi - x_mean)**2.
  !!#!!!#!! = = = = = = = = = = = = = = = = = = = = = = = = = =

  !!#!file_loop2: do nn = 1, i_file_length

  !!#!  ! Open the current colocated file
  !!#!  read(io8, '(A)', iostat = istatus) single_file_date
  !!#!  if(istatus > 0) then
  !!#!    write(*,*) 'ERROR: Problem reading line from '//trim(file_name_file)
  !!#!    cycle file_loop2
  !!#!  else
  !!#!    write(*,*) nn, single_file_date

  !!#!    total_file_name = '/home/blake.sorenson/OMI/arctic_comp/comp_data/'&
  !!#!        //'colocated_subset_'//single_file_date//'.hdf5'

  !!#!    ! Read data from this file
  !!#!    ! ------------------------
  !!#!    call h5fopen_f(total_file_name, H5F_ACC_RDWR_F, file_id, error)
  !!#!    if(error /= 0) then
  !!#!      write(*,*) 'FATAL ERROR: Could not open file '//total_file_name
  !!#!      cycle file_loop2
  !!#!    endif

  !!#!    call read_coloc_MODIS_CLD(file_id) 
  !!#!    call read_coloc_MODIS_CH7(file_id) 
  !!#!    call read_coloc_CERES_SWF(file_id) 
  !!#!    call read_coloc_NSIDC_ICE(file_id) 
  !!#!    call read_coloc_OMI_AI(file_id) 
  !!#!    call read_coloc_OMI_SZA(file_id) 

  !!#!    !!#!minii = 9
  !!#!    !!#!maxii = 11
  !!#!    !!#!minjj = 1
  !!#!    !!#!maxjj = 2
  !!#!    !!#!write(*,*) minii, maxii, minjj, maxjj

  !!#!    !!#!write(*,*) OMI_AI_data(minjj:maxjj, minii:maxii)
  !!#!    !!#!write(*,*) sum(OMI_AI_data(minjj:maxjj, minii:maxii)), &
  !!#!    !!#!    size(OMI_AI_data(minjj:maxjj, minii:maxii)), &
  !!#!    !!#!    sum(OMI_AI_data(minjj:maxjj, minii:maxii)) /  &
  !!#!    !!#!    size(OMI_AI_data(minjj:maxjj, minii:maxii))

  !!#!    omi_loop21: do ii = 1, OMI_AI_dims(2)

  !!#!      ! Check if the current pixel is along an edge
  !!#!      if(ii == 1) then
  !!#!        minii = 1
  !!#!        maxii = ii + 1
  !!#!      else if(ii == OMI_AI_dims(2)) then
  !!#!        minii = ii - 1
  !!#!        maxii = ii
  !!#!      else
  !!#!        minii = ii - 1
  !!#!        maxii = ii + 1
  !!#!      endif

  !!#!      omi_loop22: do jj = 1, OMI_AI_dims(1)

  !!#!        write(*,*) OMI_AI_data(jj,ii)

  !!#!        ! Make sure missing AI values are not included
  !!#!        ! --------------------------------------------
  !!#!        if( (abs(OMI_AI_data(jj,ii) ) < 20) .and. &
  !!#!            (abs(CERES_SWF_data(jj,ii) ) < 1200) ) then

  !!#!          ! Identify the SZA and ICE values
  !!#!          ! -------------------------------
  !!#!          sza_idx = minloc(abs(sza_bins - OMI_SZA_data(jj,ii)), dim = 1) 
  !!#!          ice_idx = minloc(abs(ice_bins - NSIDC_data(jj,ii)), dim = 1) 

  !!#!          ! Apply the cloud mask check to determine if this pixel is cloudy or 
  !!#!          ! clear
  !!#!          !
  !!#!          if(MODIS_CLD_data(jj,ii) == 3.) then
  !!#!            ! set cloud_index to 0 and continue
  !!#!            cld_idx = 1    
  !!#!            
  !!#!          else
  !!#!            ! Calculate the average of the AI values around this pixel
  !!#!            ! --------------------------------------------------------
  !!#!            ! Check if the current pixel is along an edge
  !!#!            if(jj == 1) then
  !!#!              minjj = 1
  !!#!              maxjj = jj + 1
  !!#!            else if(jj == OMI_AI_dims(1)) then
  !!#!              minjj = jj - 1
  !!#!              maxjj = jj
  !!#!            else
  !!#!              minjj = jj - 1
  !!#!              maxjj = jj + 1
  !!#!            endif


  !!#!            avg_ai = sum(OMI_AI_data(minjj:maxjj, minii:maxii)) /  &
  !!#!                    size(OMI_AI_data(minjj:maxjj, minii:maxii))
  !!#!          
  !!#!            ! Check the AI threshold to see if this is a significantly smoky 
  !!#!            ! region
  !!#!            if( (OMI_AI_data(jj,ii) > 1.5 ) .and. ( avg_ai > 1.5 )) then

  !!#!              ! Use the NSIDC_ICE value to figure out which threshold
  !!#!              ! value to use
  !!#!              ! -----------------------------------------------------
  !!#!              if(NSIDC_data(jj,ii) <= 20.) then
  !!#!                local_thresh = ch7_thresh(1)
  !!#!              else if((NSIDC_data(jj,ii) > 20.) .and. &
  !!#!                      (NSIDC_data(jj,ii) < 80)) then
  !!#!                local_thresh = ch7_thresh(4)
  !!#!              else if((NSIDC_data(jj,ii) >= 80.) .and. &
  !!#!                      (NSIDC_data(jj,ii) <= 100)) then
  !!#!                local_thresh = ch7_thresh(2)
  !!#!              else if(NSIDC_data(jj,ii) == 254.) then
  !!#!                local_thresh = ch7_thresh(3)
  !!#!              else
  !!#!                write(*,*) 'INVALID NSIDC', NSIDC_data(jj,ii)
  !!#!              endif
  !!#!               

  !!#!              ! Check the MODIS 2.1 μm reflectance to see if this is a 
  !!#!              ! misclassification
  !!#!              if((MODIS_CH7_data(jj,ii) < local_thresh) .and. &
  !!#!                 (MODIS_CH7_data(jj,ii) /= -999.)) then
  !!#!                ! MISCLASSIFICATION of cloudy as clear. Set cloud_index to 
  !!#!                ! 0 and continue
  !!#!                cld_idx = 1
  !!#!                !write(*,*) 'MISCLASSIFICATION', MODIS_CH7_data(jj,ii), &
  !!#!                !  OMI_AI_data(jj,ii), MODIS_CLD_data(jj,ii)
  !!#!              
  !!#!              else
  !!#!                ! Correct class. of cloudy as cloudy. Set cloud_index to 1 
  !!#!                ! and continue
  !!#!                cld_idx = 2
  !!#!                !write(*,*) 'GOOD CLASSIFICATN', MODIS_CH7_data(jj,ii), &
  !!#!                !  OMI_AI_data(jj,ii), MODIS_CLD_data(jj,ii)
  !!#!           
  !!#!              endif    
  !!#!          
  !!#!            else 
  !!#!              ! The smoke here is not very high, so trust the initial cloud mask
  !!#!              ! and set cloud index to 'cloudy' 1
  !!#!              cld_idx = 2    
  !!#!              !write(*,*) 'GOOD CLASSIFICATN', MODIS_CH7_data(jj,ii), &
  !!#!              !  OMI_AI_data(jj,ii), MODIS_CLD_data(jj,ii)
  !!#!          
  !!#!            endif
  !!#!          endif

  !!#!          ! After determining which bins this pixel falls into, update the
  !!#!          ! regression variables
  !!#!          ! --------------------------------------------------------------
  !!#!          y_hat = regress_slopes(cld_idx, sza_idx, ice_idx) * OMI_AI_data(jj,ii) + &
  !!#!              regress_intercepts(cld_idx, sza_idx, ice_idx)
  !!#!          sum_pred_err_square(cld_idx, sza_idx, ice_idx) = &
  !!#!              sum_pred_err_square(cld_idx, sza_idx, ice_idx) + &
  !!#!              (CERES_SWF_data(jj,ii) - y_hat)**2. 
  !!#!          sum_x_pert_square(cld_idx, sza_idx, ice_idx) = &
  !!#!              sum_x_pert_square(cld_idx, sza_idx, ice_idx) + &
  !!#!              (OMI_AI_data(jj,ii) - meanx_values(cld_idx, sza_idx, ice_idx))**2. 


  !!#!        endif
  !!#!      enddo omi_loop22
  !!#!    enddo omi_loop21

  !!#!    ! Deallocate all arrays for the next pass
  !!#!    ! ---------------------------------------
  !!#!    call clear_arrays

  !!#!    ! Close file
  !!#!    ! ----------
  !!#!    call h5fclose_f(file_id, error)
  !!#!    if(error /= 0) then
  !!#!      write(*,*) 'FATAL ERROR: Could not close file'
  !!#!      return
  !!#!    endif

  !!#!  endif
  !!#!enddo file_loop2

  !!#!! = = = = = = = = = = = = = = = = = = = = = = = = = =
  !!#!!
  !!#!! Step 5: For each of these conditions, use the accumulated parameters to
  !!#!!         calculate the slope standard error
  !!#!!         
  !!#!! = = = = = = = = = = = = = = = = = = = = = = = = = =
  !!#!!regress_slope_error = ( (1. / (count_value - 2.)) * &
  !!#!!  (sum_pred_err_square / sum_x_pert_square)  ) ** 0.5  

  !!#!write(*,*)
  !!#!write(*,*) ' = = = = = = = = = = = = = = = = = = = = = = = = = = = '
  !!#!write(*,*)
  !!#!write(*,*) ' CALCULATING SLOPE ERROR'
  !!#!write(*,*)
  !!#!write(*,*) ' = = = = = = = = = = = = = = = = = = = = = = = = = = = '

  !!#!do ii = 1, len_ice
  !!#!  do jj = 1, len_sza
  !!#!    do nn = 1, len_cld

  !!#!      if((count_values(nn,jj,ii) /= 0) .and. &
  !!#!         (count_values(nn,jj,ii) /= 2)) then

  !!#!        regress_slope_error(nn,jj,ii) = ( ( 1. / &
  !!#!              (count_values(nn,jj,ii) - 2.)) * &
  !!#!              (sum_pred_err_square(nn,jj,ii) / &
  !!#!               sum_x_pert_square(nn,jj,ii)) ) ** 0.5  

  !!#!        write(*,*) nn,jj,ii, count_values(nn,jj,ii), &
  !!#!          regress_slope_error(nn,jj,ii)

  !!#!      endif

  !!#!    enddo
  !!#!  enddo
  !!#!enddo

  ! Close the file name file
  close(io8)


  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !
  ! Write values to an output file
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  out_file_name = 'force_effic_values.h5'

  call write_force_output_file(out_file_name)

  ! Open output file
  ! ----------------
  call h5fcreate_f(trim(out_file_name), H5F_ACC_TRUNC_F, out_file_id, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open output file'
    write(*,*) '             Not writing to file'
  else
   
    !!#!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
    !!#!!
    !!#!! Write regression slopes
    !!#!!
    !!#!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
    !!#!rank = 4
    !!#!regress_dims = (/len_cod, len_cld, len_sza, len_ice/)
    !!#!!rank = 3
    !!#!!regress_dims = (/len_cld, len_sza, len_ice/)

    !!#!! Create the dataspace
    !!#!! --------------------
    !!#!call h5screate_simple_f(rank, regress_dims, dspace_id, error)
    !!#!if(error /= 0) then
    !!#!  write(*,*) 'FATAL ERROR: could not open dataspace'
    !!#!  return
    !!#!endif

    !!#!! Create the dataset
    !!#!! ------------------
    !!#!call h5dcreate_f(out_file_id, 'regress_slopes', H5T_NATIVE_REAL, &
    !!#!                  dspace_id, dset_id, error)
    !!#!if(error /= 0) then
    !!#!  write(*,*) 'FATAL ERROR: could not open dataset'
    !!#!  return
    !!#!endif

    !!#!! Write to the dataset
    !!#!! --------------------
    !!#!call h5dwrite_f(dset_id, H5T_NATIVE_REAL, regress_slopes, &
    !!#!                 regress_dims, error)
    !!#!if(error /= 0) then
    !!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
    !!#!  return
    !!#!endif

    !!#!! Close the dataset
    !!#!! -------------------------------
    !!#!call h5dclose_f(dset_id, error)

    !!#!! Close access to data space rank
    !!#!! -------------------------------
    !!#!call h5sclose_f(dspace_id, error)

    !!#!write(*,*) 'Wrote regress_slopes'

    !!#!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
    !!#!!
    !!#!! Write regression intercepts
    !!#!!
    !!#!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  

    !!#!! Create the dataspace
    !!#!! --------------------
    !!#!call h5screate_simple_f(rank, regress_dims, dspace_id, error)
    !!#!if(error /= 0) then
    !!#!  write(*,*) 'FATAL ERROR: could not open dataspace'
    !!#!  return
    !!#!endif

    !!#!! Create the dataset
    !!#!! ------------------
    !!#!call h5dcreate_f(out_file_id, 'regress_intercepts', H5T_NATIVE_REAL, &
    !!#!                  dspace_id, dset_id, error)
    !!#!if(error /= 0) then
    !!#!  write(*,*) 'FATAL ERROR: could not open dataset'
    !!#!  return
    !!#!endif

    !!#!! Write to the dataset
    !!#!! --------------------
    !!#!call h5dwrite_f(dset_id, H5T_NATIVE_REAL, regress_intercepts, &
    !!#!                 regress_dims, error)
    !!#!if(error /= 0) then
    !!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
    !!#!  return
    !!#!endif

    !!#!! Close the dataset
    !!#!! -------------------------------
    !!#!call h5dclose_f(dset_id, error)

    !!#!! Close access to data space rank
    !!#!! -------------------------------
    !!#!call h5sclose_f(dspace_id, error)

    !!#!write(*,*) 'Wrote regress_intercepts'


    !!#!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
    !!#!!
    !!#!! Write regression slope errors
    !!#!!
    !!#!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  

    !!#!! Create the dataspace
    !!#!! --------------------
    !!#!call h5screate_simple_f(rank, regress_dims, dspace_id, error)
    !!#!if(error /= 0) then
    !!#!  write(*,*) 'FATAL ERROR: could not open dataspace'
    !!#!  return
    !!#!endif

    !!#!! Create the dataset
    !!#!! ------------------
    !!#!call h5dcreate_f(out_file_id, 'regress_slope_error', H5T_NATIVE_REAL, &
    !!#!                  dspace_id, dset_id, error)
    !!#!if(error /= 0) then
    !!#!  write(*,*) 'FATAL ERROR: could not open dataset'
    !!#!  return
    !!#!endif

    !!#!! Write to the dataset
    !!#!! --------------------
    !!#!call h5dwrite_f(dset_id, H5T_NATIVE_REAL, regress_slope_error, &
    !!#!                 regress_dims, error)
    !!#!if(error /= 0) then
    !!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
    !!#!  return
    !!#!endif

    !!#!! Close the dataset
    !!#!! -------------------------------
    !!#!call h5dclose_f(dset_id, error)

    !!#!! Close access to data space rank
    !!#!! -------------------------------
    !!#!call h5sclose_f(dspace_id, error)

    !!#!write(*,*) 'Wrote regress_slope_errors'

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
    !
    ! Write regression r2 coefficients
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  

    ! Create the dataspace
    ! --------------------
    call h5screate_simple_f(rank, regress_dims, dspace_id, error)
    if(error /= 0) then
      write(*,*) 'FATAL ERROR: could not open dataspace'
      return
    endif

    ! Create the dataset
    ! ------------------
    call h5dcreate_f(out_file_id, 'regress_r2', H5T_NATIVE_REAL, &
                      dspace_id, dset_id, error)
    if(error /= 0) then
      write(*,*) 'FATAL ERROR: could not open dataset'
      return
    endif

    ! Write to the dataset
    ! --------------------
    call h5dwrite_f(dset_id, H5T_NATIVE_REAL, regress_r2, &
                     regress_dims, error)
    if(error /= 0) then
      write(*,*) 'FATAL ERROR: could not write to dataset'
      return
    endif

    ! Close the dataset
    ! -------------------------------
    call h5dclose_f(dset_id, error)

    ! Close access to data space rank
    ! -------------------------------
    call h5sclose_f(dspace_id, error)

    write(*,*) 'Wrote regress_r2'

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
    !
    ! Write bin counts
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  

    ! Create the dataspace
    ! --------------------
    call h5screate_simple_f(rank, regress_dims, dspace_id, error)
    if(error /= 0) then
      write(*,*) 'FATAL ERROR: could not open dataspace'
      return
    endif

    ! Create the dataset
    ! ------------------
    call h5dcreate_f(out_file_id, 'regress_counts', H5T_NATIVE_REAL, &
                      dspace_id, dset_id, error)
    if(error /= 0) then
      write(*,*) 'FATAL ERROR: could not open dataset'
      return
    endif

    ! Write to the dataset
    ! --------------------
    call h5dwrite_f(dset_id, H5T_NATIVE_REAL, count_values, &
                     regress_dims, error)
    if(error /= 0) then
      write(*,*) 'FATAL ERROR: could not write to dataset'
      return
    endif

    ! Close the dataset
    ! -------------------------------
    call h5dclose_f(dset_id, error)

    ! Close access to data space rank
    ! -------------------------------
    call h5sclose_f(dspace_id, error)

    write(*,*) 'Wrote regress_counts'


    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
    !
    !  NEW VARIABLE
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  

    ! Create the dataspace
    ! --------------------
    call h5screate_simple_f(rank, regress_dims, dspace_id, error)
    if(error /= 0) then
      write(*,*) 'FATAL ERROR: could not open dataspace'
      return
    endif

    ! Create the dataset
    ! ------------------
    call h5dcreate_f(out_file_id, 'regress_counts', H5T_NATIVE_REAL, &
                      dspace_id, dset_id, error)
    if(error /= 0) then
      write(*,*) 'FATAL ERROR: could not open dataset'
      return
    endif

    ! Write to the dataset
    ! --------------------
    call h5dwrite_f(dset_id, H5T_NATIVE_REAL, count_values, &
                     regress_dims, error)
    if(error /= 0) then
      write(*,*) 'FATAL ERROR: could not write to dataset'
      return
    endif

    ! Close the dataset
    ! -------------------------------
    call h5dclose_f(dset_id, error)

    ! Close access to data space rank
    ! -------------------------------
    call h5sclose_f(dspace_id, error)

    write(*,*) 'Wrote regress_counts'


  !!#!allocate(sza_bins(len_sza))
  !!#!allocate(ice_bins(len_ice))
  !!#!allocate(cod_bins(len_cod))
  !!#!!allocate(cld_bins(len_ice))
  !!#!!allocate(cld_bins(len_cld))

  !!#!allocate(sza_edges(len_sza + 1))
  !!#!allocate(ice_edges(len_ice + 1))
  !!#!allocate(cod_edges(len_cod + 1))



    ! Close the output file
    ! ---------------------
    call h5fclose_f(out_file_id, error)
    if(error /= 0) then
      write(*,*) 'FATAL ERROR: could not close output file'
      return
    endif
  endif

  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !
  ! Deallocate arrays and close program
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  call clear_force_eff_arrays
 
  deallocate(sza_bins)
  deallocate(ice_bins)
  deallocate(cod_bins)
  deallocate(sza_edges)
  deallocate(ice_edges)
  deallocate(cod_edges)

  ! Close the HDF5 interface
  ! ------------------------
  call h5close_f(error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: Could not close HDF5 library'
    return
  endif
  write(*,*) 'Interface closed'

end program test_force_eff_calc
