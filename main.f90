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

    REAL*8            a(2,Lmax+1,Lmax+1)    !Coefficients
    REAL*8            S(Nth,Nphi)           !Temperatures

    ! Parameters for FFT using Gauss-Legendre Quadrature
    REAL*8            GLQWeights(Lgrid+1)  !Gauss-Legendre Quadrature weights
    REAL*8            GLQZeros(Lgrid+1)    !Gauss-Legendre Quadrature zeros
    
    REAL*8            walltime_start, walltime_stop  !calc. total wall time
    REAL*8            cputime_start, cputime_stop    !calc. program speed
    REAL*8            eps         !Machine epsilon (needs to be pre-computed)

    REAL*8              R               !Radius of the object
    REAL*8              theta(Nth)      !The angle theta
    REAL*8              phi(Nphi)       !The angle phi

    REAL*8              T0          !Initial temperature
    REAL*8              tolA        !Tolerance of the solution
    REAL*8              tolF        !Tolerance of the residual

    INTEGER*4           l, ml       !Degree of spherical harmonics
    INTEGER*4           LWORK       !For inverting Jacobian
 
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
    R           = 0.2D0                 !Radius of the object in meters
    !Tolerances
    tolA        = 1000.0D0              !Tolerance of solution
    tolF        = 0.1D0                 !Tolerance of residual

    LWORK       = 18496

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

    CALL GetInitialData(&
    & Nth, Nphi,&
    & Lmax, Lgrid,&
    & GLQWeights, GLQZeros,&
    & T0,&
    & R, theta, phi,&
    & a)

    CALL SolveEquation(&
    & Nth, Nphi, Lmax, Lgrid, GLQWeights, GLQZeros,&
    & R, theta, phi,&
    & Maxit, tolA, tolF,&
    & LWORK,&
    & a) 

    CALL GetResults(&
    & Nth, Nphi,&
    & Lmax, Lgrid,&
    & GLQWeights, GLQZeros,&
    & R, theta, phi,&
    & a, S)

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


