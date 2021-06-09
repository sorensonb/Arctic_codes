program omi_frequency
!
! NAME:
!   omi_frequency.f90
!
! PURPOSE:
! 
! CALLS:
!   mie_calc.f90
!
! MODIFICATIONS:
!   Blake Sorenson <blake.sorenson@und.edu>     - 2018/10/24:
!     Written
!
!  ############################################################################

  implicit none

  complex                :: m            ! refractive index (real and imaginary
                                         ! components.
  real                   :: diam         ! Particle diameter (used to find size
                                         ! parameter) (microns)
  real                   :: lambda       ! wavelength  (microns)
  real,dimension(3)      :: reals        ! real components of refractive index
  real,dimension(3)      :: size_params  ! size parameters
  real                   :: pi
 
  integer                :: i            ! loop counter

  ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  ! Set up count variables to count the number of grid boxes with
  ! high AI values

  ! Loop over file names, once hour of file name goes outside the
  ! +/- 3 hrs from synoptic time, print counts and date

   
  
  

  pi = 3.14159265
  lambda = 3.7 ! in microns

  ! Put the real components of the refractive indices used in Table 1 of
  ! Wang and van de Hulst (1991) into an array which will be accessed
  ! in the following loops.
  reals = [1.55,1.50,1.342] 
  ! These size parameters are the ones used in the table.
  size_params = [5.21282,100.,1570.7963]

  write(*,*)
  write(*,*) 'Non-absorbing refractive index m:'
  write(*,*) '________________________________________________________________'

  do i=1,3
    m = complex(reals(i),0.00)
    diam   = (size_params(i)*lambda)/pi
    write(*,*) 'diameter (microns) =',diam
    write(*,*) 'wavelength (microns) =',lambda
    call mie_calc(m,lambda,diam)

  enddo

  write(*,*)
  write(*,*) 'Absorbing refractive index m:'
  write(*,*) '________________________________________________________________'

  do i=1,3
    m = complex(reals(i),-0.1)
    diam   = (size_params(i)*lambda)/pi
    write(*,*) 'diameter (microns) =',diam
    write(*,*) 'wavelength (microns) =',lambda
    call mie_calc(m,lambda,diam)

  enddo

end program omi_frequency
