subroutine mie_calc(m,lambda,diam)
!
!  NAME:
!    mie_calc.f90
!
!  PURPOSE:
!    Calculate Q_ext, Q_sca, Q_abs, and g 
!
!  CALLS:
!    bessel.f90, newmann.f90, calc_deriv_A.f90, calc_deriv_B.f90
!
!  MODIFICATIONS:
!    Blake Sorenson <blake.sorenson@und.edu>    - 2018/10/23:
!      Written
!
!  ###########################################################################

  implicit none

  real                             :: lambda       ! wavelength 
  real                             :: diam         ! particle diameter
  real                             :: m_real       ! real component of
                                                   ! refractive index
  real                             :: m_imag       ! imaginary component of
                                                   ! refractive index

  real                             :: pi           ! pi

  real                             :: Q_sca        ! scattering coefficient
  real                             :: Q_ext        ! extinction coefficient
  real                             :: Q_abs        ! absorption coefficient
  real                             :: g_asym       ! asymmetry factor
  real                             :: g_asym_p2    ! part 2 of the asymmetry
                                                   ! factor calculation.

  integer                          :: i             ! loop counter
  integer                          :: n_max         ! the maximum number of
                                                    ! iterations for Mie 
                                                    ! series convergence

  complex                          :: rho           ! size parameter (used
                                                    ! in Mie calculations)
  complex                          :: eta           ! complex argument used
                                                    ! in mie calculations
                                                    ! = m*rho

  complex                          :: psi           ! Bessel function value
  complex                          :: chi           ! Newmann function value
  complex                          :: imaginary     ! Imaginary value (i)
  complex                          :: m             ! Refractive index
  complex                          :: k             ! 
  complex                          :: a_coeff       ! Mie coefficient a
  complex                          :: b_coeff       ! Mie coefficient b
  complex                          :: deriv_A       ! derivative A with eta
                                                    ! (psi_prime/psi)                                                     
  complex                          :: deriv_B       ! derivative B
                                                    ! (chi_prime/chi)
  complex                          :: deriv_A_rho   ! derivative A with rho
                                                    
  complex,dimension(:),allocatable :: an_array      ! array to hold the a
                                                    ! coefficients
  complex,dimension(:),allocatable :: bn_array      ! array to hold the b 
                                                    ! coefficients

  ! Declare function types
  complex                          :: bessel
  complex                          :: newmann
  complex                          :: calc_deriv_B
  complex                          :: calc_deriv_A

  ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
   
  ! Define the refractive index
  pi = 3.14159265
  imaginary = (0,1)
 
  ! Set up the size parameter 
  k = complex(pi*diam/lambda,0) ! x, !(pi*diam)/lambda 

  ! Set up the variables to use in the Mie calculations
  rho = k
  eta = k*m

  ! Calculate the maximum number of iterations for the Mie code to converge
  n_max = int(rho+4.*(rho**(1.0/3.0))+2.)

  ! Allocate space for the an and  bn terms 
  allocate(an_array(n_max))
  allocate(bn_array(n_max))

  ! Set the initial values of these values to 0 for summation.
  Q_sca = 0.
  Q_ext = 0.
  g_asym = 0.
  g_asym_p2 = 0.

  ! Calculate the Mie coefficients
  do i=1,n_max
    ! Call subroutines to calculate the values of the Bessel and Newmann
    ! functions.
    psi = bessel(rho,real(i))
    chi = newmann(rho,real(i))

    ! Calculate the derivative values that will be used to calculate the
    ! an and bn values.
    deriv_A = calc_deriv_A(eta,real(i))
    deriv_B = calc_deriv_B(rho,real(i))
    deriv_A_rho = calc_deriv_A(rho,real(i))

    ! Calculate the an and bn Mie coefficients.
    a_coeff = 1./(1.+(imaginary*(chi*(deriv_A-m*deriv_B)))/                   &
    (psi*(deriv_A-m*deriv_A_rho)))
    b_coeff = 1./(1.+(imaginary*(chi*(m*deriv_A-deriv_B)))/                   &
    (psi*(m*deriv_A-deriv_A_rho)))

    ! Put the Mie coefficients into arrays to be used later to find the
    ! asymmetry factor.
    an_array(i) = a_coeff
    bn_array(i) = b_coeff

    ! Calculate the summations for Q_sca and Q_ext
    Q_sca = Q_sca+(2.*i+1.)*(abs(a_coeff)**2.+abs(b_coeff)**2.)
    Q_ext = Q_ext+(2.*i+1.)*real(a_coeff+b_coeff)
    g_asym_p2 = g_asym_p2+((2.*i+1.)/(i*(i+1.)))*real(a_coeff*conjg(b_coeff))
  enddo

  ! Calculate Q_sca, Q_ext, and g
  Q_sca = 2*Q_sca/(real(k)**2.)
  Q_ext = 2*Q_ext/(real(k)**2.)
  Q_abs = Q_ext-Q_sca

  ! Calculate the asymmetry factor
  do i=1,n_max-1
    g_asym = g_asym+((i*(i+2.))/(i+1.))*real((an_array(i)*                    &
             conjg(an_array(i+1)))+(bn_array(i)*conjg(bn_array(i+1))))
  enddo
  g_asym = (4./(real(k)**2.*Q_sca))*(g_asym+g_asym_p2)

  write(*,*) 'x = ',real(k),'                       Q_ext   =',Q_ext
  write(*,*) 'm = ',m,      '     Q_abs   =',Q_abs
  write(*,*) '                                            Q_sca   =',Q_sca
  write(*,*) '                                            g_asym  =',g_asym

  ! Free allocated arrays
  deallocate(an_array)
  deallocate(bn_array)

end subroutine mie_calc
