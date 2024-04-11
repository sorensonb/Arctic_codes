module force_eff_vars
!
! NAME:
!   force_eff_vars
!
! PURPOSE:
!   Contain the data arrays and dimensions for all variable types, as well as
!   some functions.
! 
! CALLS:
!   None.
!
! MODIFICATIONS:
!   Blake Sorenson <blake.sorenson@und.edu>     - 2021/07/09:
!     Written
!
!  ############################################################################

  use hdf5

  implicit none

  integer               :: len_sza
  integer               :: len_ice
  integer               :: len_cld
  integer               :: len_cod

  real, dimension(:), allocatable           :: sza_bins
  real, dimension(:), allocatable           :: ice_bins
  real, dimension(:), allocatable           :: cod_bins

  real, dimension(:), allocatable           :: sza_edges
  real, dimension(:), allocatable           :: ice_edges
  real, dimension(:), allocatable           :: cod_edges


  real, dimension(:,:,:,:), allocatable     :: sumx_values
  real, dimension(:,:,:,:), allocatable     :: sumy_values
  real, dimension(:,:,:,:), allocatable     :: sumx2_values
  real, dimension(:,:,:,:), allocatable     :: sumy2_values
  real, dimension(:,:,:,:), allocatable     :: sumxy_values
  integer, dimension(:,:,:,:), allocatable  :: count_values

  real, dimension(:,:,:,:), allocatable     :: regress_slopes
  real, dimension(:,:,:,:), allocatable     :: regress_intercepts
  real, dimension(:,:,:,:), allocatable     :: regress_slope_error
  real, dimension(:,:,:,:), allocatable     :: regress_r2


  contains
 
    subroutine allocate_force_out_arrays(dim1, dim2, dim3, dim4)
    ! --------------------------------------------------------------
    ! This subroutine allocates the output data and dimension arrays
    ! --------------------------------------------------------------

    integer   :: dim1
    integer   :: dim2
    integer   :: dim3
    integer   :: dim4
    integer   :: error

    ! Dimensions for output data
    ! --------------------------

    allocate(sumx_values(dim1, dim2, dim3, dim4), stat = error)
    allocate(sumy_values(dim1, dim2, dim3, dim4), stat = error)
    allocate(sumx2_values(dim1, dim2, dim3, dim4), stat = error)
    allocate(sumy2_values(dim1, dim2, dim3, dim4), stat = error)
    allocate(sumxy_values(dim1, dim2, dim3, dim4), stat = error)
    allocate(count_values(dim1, dim2, dim3, dim4), stat = error)
    allocate(regress_slopes(dim1, dim2, dim3, dim4), stat = error)
    allocate(regress_intercepts(dim1, dim2, dim3, dim4), stat = error)
    allocate(regress_slope_error(dim1, dim2, dim3, dim4), stat = error)
    allocate(regress_r2(dim1, dim2, dim3, dim4), stat = error)

    sumx_values(:,:,:,:) = 0.
    sumy_values(:,:,:,:) = 0.
    sumx2_values(:,:,:,:) = 0.
    sumy2_values(:,:,:,:) = 0.
    sumxy_values(:,:,:,:) = 0.
    count_values(:,:,:,:) = 0
    regress_slopes(:,:,:,:) = 0.
    regress_intercepts(:,:,:,:) = 0.
    regress_slope_error(:,:,:,:) = 0.
    regress_r2(:,:,:,:) = 0.


    if ( error < 0 ) then
       write(*,*) " *** Error allocating H5datasets in output data"
       return
    endif
      
    write(*,*) "Allocated out data to dimensions",dim1,dim2,dim3,dim4

    end subroutine allocate_force_out_arrays
 
    subroutine clear_force_eff_arrays
    ! -------------------------------------------------------------
    ! This subroutine deallocates all the data and dimension arrays
    ! -------------------------------------------------------------
    
      ! ############################################################################


      ! Dimensions for input data
      ! ------------------------- 

      if(allocated(regress_slopes))       deallocate(regress_slopes)
      if(allocated(regress_intercepts))   deallocate(regress_intercepts)
      if(allocated(regress_slope_error))  deallocate(regress_slope_error)
      if(allocated(regress_r2))           deallocate(regress_r2)
                                        
      if(allocated(sumx_values))          deallocate(sumx_values)
      if(allocated(sumy_values))          deallocate(sumy_values)
      if(allocated(sumx2_values))         deallocate(sumx2_values)
      if(allocated(sumy2_values))         deallocate(sumy2_values)
      if(allocated(sumxy_values))         deallocate(sumxy_values)
      if(allocated(count_values))         deallocate(count_values)
                                        
      if(allocated(sza_bins))             deallocate(sza_bins)
      if(allocated(ice_bins))             deallocate(ice_bins)
      if(allocated(cod_bins))             deallocate(cod_bins)
      if(allocated(sza_edges))            deallocate(sza_edges)
      if(allocated(ice_edges))            deallocate(ice_edges)
      if(allocated(cod_edges))            deallocate(cod_edges)

    end subroutine clear_force_eff_arrays

end module force_eff_vars
