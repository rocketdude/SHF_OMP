!======================================================! 
! Radiation Boundary Condition Solver                  ! 
!======================================================!

  PROGRAM              RadiationBoundary
    USE                omp_lib
    USE                SHTOOLS
    USE                HDF5
    IMPLICIT           none

    INTEGER*4, PARAMETER ::        Lmax     = 16
    INTEGER*4, PARAMETER ::        Lgrid    = 16

    INTEGER*4, PARAMETER ::        Nth      = Lgrid+1
    INTEGER*4, PARAMETER ::        Nphi     = 2*Lgrid+1

    REAL*8, PARAMETER ::    PI = 3.141592653589793238462643383279502884197D0

    !Lmax is the maximum value of l considered
    !Lgrid is the maximum value of l that can be found from 

!-------------------------------------------------------!
!     Declare Commons                                   !
!-------------------------------------------------------!

!-------------------------------------------------------!
!     DATA I/O                                          !
!-------------------------------------------------------!

    CHARACTER*32      CTemp
    LOGICAL           FileExist

!--------------------------------------------------------!
!     Declare Numerical Inputs                           ! 
!--------------------------------------------------------!

    INTEGER*4           it     !Iteration counter
    INTEGER*4           Maxit  !Maximum iteration

!--------------------------------------------------------!
!     Declare Level Set Data                             !
!--------------------------------------------------------!

    ! Things that need to be allocated due to sph. harmonic FFTW
    COMPLEX*16        a(2,Lmax+1,Lmax+1) 

    ! Parameters for FFT using Gauss-Legendre Quadrature
    REAL*8            GLQWeights(Lgrid+1)  !Gauss-Legendre Quadrature weights
    REAL*8            GLQZeros(Lgrid+1)    !Gauss-Legendre Quadrature zeros
    
    REAL*8            walltime_start, walltime_stop  !calc. total wall time
    REAL*8            cputime_start, cputime_stop    !calc. program speed
    REAL*8            eps         !Machine epsilon (needs to be pre-computed)

    REAL*8              theta(Nth)      !The angle theta
    REAL*8              phi(Nphi)        !The angle phi

    REAL*8              T0          !Initial temperature
    REAL*8              SolPhi      !Solar constant
    REAL*8              epsStar     !Normalized emissivity
    REAL*8              alpStar     !Normalized absorptivity
    REAL*8              kappaStar   !Normalized conductivity

    INTEGER*4           l, ml       !Degree of spherical harmonics
 
!--------------------------------------------------------!
!     Declare Local Parameters                           !
!--------------------------------------------------------!

    INTEGER*4           i,j,k       !Counters for DO loops
    INTEGER*4           nthreads    !Number of threads

!========================================================!
!     MAIN PROGRAM                                       !
!========================================================!

!--------------------------------------------------------!
!     Parameters                                         !
!--------------------------------------------------------!
    
    !Parameters related to initial conditions

    !Termination conditions
    Maxit = 1000

    !Parameters related to the heat problem
    T0          = 300.0D0               !Initial temperature guess in K
    SolPhi      = 1366.0D0              !Solar constant in W/m^2
    epsStar     = 0.90D0/2.0D2          !Normalized emissivity (eps/kappa)
    alpStar     = 0.45D0/2.0D2          !Normalized absorptivity (alpha/kappa)
    kappaStar   = 2.0D2/(1.5D4*2.0D2)   !kappa / rho*cp

    eps =  2.22044604925031308D-016 !Machine epsilon (precalculate)
 
!--------------------------------------------------------!
!     Timer Start                                        !
!--------------------------------------------------------!	  

    CALL cpu_time( cputime_start )
    walltime_start = OMP_GET_WTIME()
!--------------------------------------------------------!
!     Compute Spherical Harmonic Transform Parameters    !
!--------------------------------------------------------!

    IF( Lmax .GT. Lgrid ) THEN
        STOP "***ERROR*** Not enough collocation points"
    END IF

    CALL PreCompute(Lgrid, GLQZeros, GLQWeights, NORM=1, CNORM=1)

!--------------------------------------------------------!
!     Creating mesh & Inquiring Threads                  !
!--------------------------------------------------------!

!!$    !Set the number of MKL threads (dynamic is automatic)
!!$    CALL MKL_SET_DYNAMIC(.TRUE.)

    !Inquire number of threads
    nthreads = omp_get_max_threads()
    WRITE(*,*) 'No. of threads = ', nthreads

    CALL GetThetaAndPhi(Nth,Nphi,Lgrid,theta,phi)
     
!--------------------------------------------------------!
!     Echo certain parameters                            !
!--------------------------------------------------------!

    PRINT *, 'epsilon star = ', epsStar
    PRINT *, 'alpha star = ', alpStar
    PRINT *, 'kappa star = ', kappaStar
    PRINT *, 'Lmax = ', Lmax
    PRINT *, '(Nth, Nphi) = ', Nth, ',', Nphi
    PRINT *, 'theta = ', theta
    PRINT *, 'phi = ', phi
    PRINT *, 'Max # of iterations =', Maxit
    
!--------------------------------------------------------!
!     Evolve Eikonal S                                   !
!--------------------------------------------------------!

    PRINT *, '==============================='
    PRINT *, 'START MAIN PROGRAM'


    PRINT *, 'END PROGRAM'
    PRINT *, '==============================='


!--------------------------------------------------------!
!     Timer Stop                                         !
!--------------------------------------------------------!

    CALL cpu_time( cputime_stop )
    walltime_stop = OMP_GET_WTIME()
    WRITE(*,*) 'Elapsed CPU Time = ', cputime_stop - cputime_start
    WRITE(*,*) 'Elapse Walltime = ', walltime_stop - walltime_start

    STOP
  END PROGRAM RadiationBoundary


