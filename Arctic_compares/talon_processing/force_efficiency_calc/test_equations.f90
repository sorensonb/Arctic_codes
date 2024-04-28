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

  implicit none


  integer               :: nn
  integer               :: num_values

  real                 :: sumx_values
  real                 :: sumy_values
  real                 :: sumx2_values
  real                 :: sumy2_values
  real                 :: sumxy_values
  integer              :: count_values

  real                 :: meanx_values
  real                 :: meany_values
  real                 :: sum_pred_err_square
  real                 :: sum_x_pert_square

  real                 :: regress_slopes
  real                 :: regress_intercepts
  real                 :: regress_slope_error
  real                 :: regress_slope_error_v2
  real                 :: regress_r2

  real                 :: s_xy
  real                 :: s_xx
  real                 :: sse 
  real                 :: sst 

  real                  :: y_hat

  real, dimension(:)    :: x_vals(10)
  real, dimension(:)    :: y_vals(10)

  ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !
  ! Step 1:  Determine how many bins for each variable are needed
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  ! Specify the bin ranges and sizes
  ! --------------------------------

  ! Allocate arrays to hold the 2.1 μm thresholds for splitting
  ! between clear and cloudy conditions
  ! The indices follow the ice type:
  !   1: ocean
  !   2: ice
  !   3: land
  !   4: mix
  ! -----------------------------------------------------------

  x_vals = [1.,2.,3.,4.,5.,6.,7.,8.,9.,10.]
  y_vals = [1.,2.,3.,4.,5.,6.,7.,8.,9.,10.]
  
  num_values = size(x_vals) 

  sumx_values  = 0.
  sumy_values  = 0.
  sumx2_values = 0.
  sumy2_values = 0.
  sumxy_values = 0.
  count_values = 0

  meanx_values = 0.
  sum_pred_err_square = 0.
  sum_x_pert_square = 0.

  regress_slopes = 0.
  regress_intercepts = 0.
  regress_slope_error = 0.
  regress_slope_error_v2 = 0.
  regress_r2 = 0.

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

  do nn = 1, num_values
    ! After determining which bins this pixel falls into, update the
    ! regression variables
    ! --------------------------------------------------------------
    sumx_values = sumx_values + x_vals(nn)
    sumy_values = sumy_values + y_vals(nn)
    sumx2_values = sumx2_values + x_vals(nn)**2.
    sumy2_values = sumy2_values + y_vals(nn)**2.
    sumxy_values = sumxy_values + (x_vals(nn) * y_vals(nn))
    count_values = count_values + 1
  enddo

  ! = = = = = = = = = = = = = = = = = = = = = = = = = =
  !
  ! Step 3: For each of these conditions, use the accumulated parameters to
  !         calculate the slope and intercept of the regression lines
  !         
  ! = = = = = = = = = = = = = = = = = = = = = = = = = =
  !!#!regress_intercepts = ((sumy_values * &
  !!#!    (sumx_values**2.)) - (sumx_values * &
  !!#!     sumxy_values)) / &
  !!#!    (count_values * sumx2_values - &
  !!#!     (sumx_values**2.))

  !!#!regress_slopes = ((count_values * &
  !!#!    sumxy_values) - (sumx_values - &
  !!#!    sumy_values)) / &
  !!#!    (count_values * sumx2_values - &
  !!#!     (sumx_values**2.))

  ! from https://www.ncl.ac.uk/webtemplate/ask-assets/external/&
  !       maths-resources/statistics/regression-and-correlation/&
  !       simple-linear-regression.html

  !!#!regress_slopes = (sumxy_values - ( (sumx_values * sumy_values) &
  !!#!   / count_values )  ) / ( sumx2_values - ( (sumx_values**2.) / count_values )  )


  ! Calculate the mean x (AI) value here per bin
  ! --------------------------------------------
  meanx_values = sumx_values / count_values
  meany_values = sumy_values / count_values

  s_xy = sumxy_values - ( (sumx_values * sumy_values) / count_values)
  s_xx = sumx2_values - ( (sumx_values**2.) / count_values ) 

  regress_slopes = s_xy / s_xx

  regress_intercepts = meany_values - regress_slopes * meanx_values

  write(*,*) 'sumx: ', sumx_values
  write(*,*) 'sumy: ', sumy_values
  write(*,*) 'sumx2:', sumx2_values
  write(*,*) 'sumxy:', sumxy_values
  write(*,*) 'count:', count_values
  write(*,*) 'meanx:', meanx_values
  write(*,*) 'meany:', meanx_values

  write(*,*) 'slope:    ', regress_slopes
  write(*,*) 'intercept:', regress_intercepts

  ! = = = = = = = = = = = = = = = = = = = = = = = = = =
  !
  ! Step 4: Loop back over all the colocated data products \
  !         and calculate, for each
  !         surface type and SZA, the following
  !         (y = ceres_swf, x = omi_ai)
  !
  !
  !         Need the following:
  !
  !       y_hat              : predicted SWF value given the regression and AI
  !       Σ(yi - y_hati)**2.
  !       Σ(xi - x_mean)**2.
  ! = = = = = = = = = = = = = = = = = = = = = = = = = =

  sse = sumy2_values - regress_intercepts * sumy_values - &
    regress_slopes * sumxy_values

  sst = sumy2_values - (sumy_values**2.) / count_values

  regress_r2 = 1. - (sse / sst)
 
  ! After determining which bins this pixel falls into, update the
  ! regression variables
  ! --------------------------------------------------------------
  !!#!do nn = 1, num_values
  !!#!  y_hat = regress_slopes * x_vals(nn) + &
  !!#!      regress_intercepts
  !!#!  sum_pred_err_square = &
  !!#!      sum_pred_err_square + &
  !!#!      (y_vals(nn) - y_hat)**2. 
  !!#!  sum_x_pert_square = &
  !!#!      sum_x_pert_square + &
  !!#!      (x_vals(nn) - meanx_values)**2. 

  !!#!  write(*,*) nn, y_vals(nn), y_hat, x_vals(nn)

  !!#!enddo

  !!#!write(*,*) 'sum_pred_err:', sum_pred_err_square
  write(*,*) 'sse         :', sse
  write(*,*) 'r2          :', regress_r2

  ! = = = = = = = = = = = = = = = = = = = = = = = = = =
  !
  ! Step 5: For each of these conditions, use the accumulated parameters to
  !         calculate the slope standard error
  !         
  ! = = = = = = = = = = = = = = = = = = = = = = = = = =
  !regress_slope_error = ( (1. / (count_value - 2.)) * &
  !  (sum_pred_err_square / sum_x_pert_square)  ) ** 0.5  

  write(*,*)
  write(*,*) ' = = = = = = = = = = = = = = = = = = = = = = = = = = = '
  write(*,*)
  write(*,*) ' CALCULATING SLOPE ERROR'
  write(*,*)
  write(*,*) ' = = = = = = = = = = = = = = = = = = = = = = = = = = = '

  regress_slope_error_v2 = ( (1. / (count_values - 2)) * ( sse / s_xx  )  ) ** 0.5

  !!#!regress_slope_error = ( ( 1. / &
  !!#!      (count_values - 2)) * &
  !!#!      (sse / &
  !!#!       sum_x_pert_square) ) ** 0.5  

  write(*,*) 'Slope error: ', regress_slope_error
  write(*,*) 'Slope error2:', regress_slope_error_v2

end program test_force_eff_calc
