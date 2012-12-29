!======================================================! 
! Spectral Horizon Finder (SHF)                        ! 
!======================================================!

! This program is an event horizon finder using the pseudospectral methods.
! The eikonal data S(t, r, theta, phi) is represented as orthogonal functions 
! of Chebyshev polynomials & spherical harmonics.
! In other words, we are approximating 
! S(t, r, theta, phi) ~= a_nlm(t) T_n (rho) Y_lm (theta, phi).
! This is an improvement to the old SHF with the use of FFT and SHTOOLS.

  PROGRAM              SHF
    USE                omp_lib
    USE                SHTOOLS
    USE                HDF5
    IMPLICIT           none

    INTEGER*4, PARAMETER ::        Mr       = 10
    INTEGER*4, PARAMETER ::        Lmax     = 2
    INTEGER*4, PARAMETER ::        Lgrid    = 4
    INTEGER*4, PARAMETER ::        TP       = 4
    INTEGER*4, PARAMETER ::        SpM      = 6

    INTEGER*4, PARAMETER ::        Nr       = Mr+1
    INTEGER*4, PARAMETER ::        Nth      = Lgrid+1
    INTEGER*4, PARAMETER ::        Nphi     = 2*Lgrid+1
    ! Mr is the degree of the radial Chebyshev polynomial
    ! Lmax is the maximum values of l and m of the spherical harmonics
    ! Nr, Nth, and Nphi are the # of points in r, theta, phi directions
    ! TP-1 is the temporal interpolation order
    ! SpM is the number of additional directions we'd like to compute U for

    REAL*8, PARAMETER ::    PI = 3.141592653589793238462643383279502884197D0

!-------------------------------------------------------!
!     Declare Commons                                   !
!-------------------------------------------------------!

!-------------------------------------------------------!
!     DATA I/O                                          !
!-------------------------------------------------------!

    CHARACTER*32      CTemp
    CHARACTER*32      aFile
    LOGICAL           FileExist

    INTEGER*4           WriteSit        
    !iteration at which we output S into file
    INTEGER*4           Writeait        
    !iteration at which we output the coef. a_nlm
    INTEGER*4           WriteRit
    !iteration at which we output the temperature along a radial line

!--------------------------------------------------------!
!     Declare Numerical Inputs                           ! 
!--------------------------------------------------------!

    INTEGER*4           it     !Iteration counter
    INTEGER*4           Startit!Starting iteration
    INTEGER*4           Maxit  !Maximum iteration

    INTEGER*4           SFLAG  
    !If SFLAG = 1, then we are continuing previous run

    REAL*8            cfl    !Courant-Friedrich-Lewy factor
    REAL*8            c      !Steepness of the initial tanh function

    INTEGER*4         STAT   !Status placeholder

!--------------------------------------------------------!
!     Declare Level Set Data                             !
!--------------------------------------------------------!

    ! We use the usual spherical coordinates
    ! X = r*sin(theta)*cos(phi)
    ! Y = r*sin(theta)*sin(phi)
    ! Z = r*cos(theta)

    ! Things that need to be allocated due to sph. harmonic FFTW
    COMPLEX*16        a(Mr+1,2,Lmax+1,Lmax+1) 
            !time-dependent coefficients a_nlm(t)

    ! Parameters for FFT using Gauss-Legendre Quadrature
    REAL*8            GLQWeights(Lgrid+1)  !Gauss-Legendre Quadrature weights
    REAL*8            GLQZeros(Lgrid+1)    !Gauss-Legendre Quadrature zeros
    
    REAL*8            t           !Current time
    REAL*8            tfinal      !Final time
    REAL*8            tdir        !Direction of time
    REAL*8            dt          !Time increment
    REAL*8            dx          !Smallest distance between two points
    REAL*8            rdtheta     !Smallest distance in the theta direction
    REAL*8            rsinthdphi  !Smallest distance in the phi direction
    REAL*8            walltime_start, walltime_stop  !calc. total wall time
    REAL*8            cputime_start, cputime_stop    !calc. program speed
    REAL*8            eps         !Machine epsilon (needs to be pre-computed)

    REAL*8            rho(Nr)   !Canonical radial coord. for Chebyshev poly.
    REAL*8            r(Nr)     !The radial coordinate r
    REAL*8            rmax, rmin  !rmax and rmin define the radial domain,
                                  !used to calculate rho

    REAL*8              theta(Nth)      !The angle theta
    REAL*8              phi(Nphi)        !The angle phi
    REAL*8              thetaSp(SpM)!The angle theta in specific directions
    REAL*8              phiSp(SpM)  !The angle phi in specific directions

    REAL*8              X0,Y0,Z0    !X0,Y0,Z0 are the initial axes of the
                                    !the spheroids
                                    !X0=Y0=Z0 is a sphere

    REAL*8              T0          !Initial temperature
    REAL*8              SolPhi      !Solar constant
    REAL*8              epsStar     !Normalized emissivity
    REAL*8              alpStar     !Normalized absorptivity

    REAL*8              TempAve     !Average temperature

    INTEGER*4           n           !Degree of Chebyshev polynomial
    INTEGER*4           l, ml       !Degree of spherical harmonics
    INTEGER*4           reinit      !Iteration at which we reinitialize
 
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
    SFLAG = 0                       !If SFLAG = 1, 
                                    !we are continuing previous run: 
                                    !change t, Startit and aFile
    t = 0.0D0                       !Last time from previous run
    Startit = 6751                  !Startit = last iteration + 1
    aFile = 'a10.dat'

    !Termination conditions
    Maxit = 5
    tfinal = 2.5D0

    IF( SFLAG .EQ. 0 ) THEN
       Startit = 1                  !Starting from iteration 1
    END IF

    X0 = 1.0D0                    !Initial axes of the spheroids in the
    Y0 = 1.0D0                    !X, Y, and Z directions
    Z0 = 1.0D0

    !Parameters related to the heat problem
    T0          = 400                   !Initial temperature in K
    SolPhi      = 1366                  !Solar constant in W/m^2
    epsStar     = 0.9/200               !Normalized emissivitity (eps/k)
    alpStar     = 0.45/200              !Normalized absorptivity (alpha/k)

    !Simulation parameters                              
    !Note: negative rootsign, positive lapse & shift functions, 
    !      and negative tdir give EH finder
    eps =  2.22044604925031308D-016 !Machine epsilon (precalculate)
    c = 0.1D0
    cfl = 1.508D0                   !Depends on which SSP-Runge-Kutta used
                                    !SSPRK(5,4)=>cfl=1.508 and 
                                    !SSPRK(3,3)=>cfl=1.0
    tdir = +1.0D0                   !Direction of time, choose +1.0D0 or -1.0D0
    reinit = 15

    !Spherical grid parameters
    rmax = 1.20D0                   !maximum value of r
    rmin = 0.20D0                   !minimum value of r

    !Additional directions we'd like to compute U
    thetaSp = (/ 0.0D0, PI, PI/2.0D0, PI/2.0D0, PI/2.0D0, PI/2.0D0 /)
    phiSp = (/ 0.0D0, 0.0D0, 0.0D0, PI/2.0D0, PI, 1.5D0*PI /)
    
!--------------------------------------------------------!
!     Output Parameters                                  !
!--------------------------------------------------------!

    WriteSit         = 100000
    Writeait         = 1000

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

    !Set up Chebyshev-Gauss-Radau points and DH collocation points
    !$OMP PARALLEL DO
    DO i = 0, Mr
       rho(i+1) = -1.0D0*COS( 2*PI*i / (2*Mr+1) ) 
       r(i+1) = 0.5D0*( (rmax + rmin) - (rmax - rmin)*rho(i+1) )
    END DO
    !$OMP END PARALLEL DO

    CALL GetThetaAndPhi(Nth,Nphi,Lgrid,theta,phi)
     
!--------------------------------------------------------!
!     Echo certain parameters                            !
!--------------------------------------------------------!

!!$    PRINT *, 'Mass of blackhole =', Mass
    PRINT *, 'Mr = ', Mr
    PRINT *, 'Lmax = ', Lmax
    PRINT *, '(Nr, Nth, Nphi) = ', Nr, ',', Nth, ',', Nphi
    PRINT *, 'theta = ', theta
    PRINT *, 'phi = ', phi
    PRINT *, '# of iterations =', (Maxit-Startit+1)
    PRINT *, 'Reinitializing every ', reinit, 'iterations'

!--------------------------------------------------------!
!     Find smallest spatial increment                    !
!--------------------------------------------------------!

    dx = rmax - rmin
    DO i = 1, Mr
       IF( (r(i+1) - r(i)) .LT. dx ) THEN
          dx = ABS(r(i+1) - r(i))
       END IF
    END DO
    
    !then set the value of dt
    dt = tdir * cfl * dx
    
!--------------------------------------------------------!
!     Initial Data                                       !
!--------------------------------------------------------!

    PRINT *, '==============================='
    PRINT *, 'INITIALIZING PROGRAM'

    IF( SFLAG .EQ. 1 ) THEN
       PRINT *, 'Continuing run'
       PRINT *, 'Previous iteration=' ,(Startit-1)

       CALL Read4dC(Mr+1, 2, Lmax+1, Lmax+1, a, aFile)

    ELSE

       CALL GetInitialData(& 
            & Nr, Nth, Nphi,&
            & Mr, Lmax, Lgrid,&
            & GLQWeights, GLQZeros,&
            & c,&
            & X0, Y0, Z0,&
            & T0,&
            & r, rho, theta, phi,&
            & a)
    END IF

    PRINT *, '==============================='


!--------------------------------------------------------!
!     Evolve Eikonal S                                   !
!--------------------------------------------------------!

    IF( SFLAG .EQ. 0 ) THEN
       !Initialize output files and write values for iteration 0
       CTemp = 'Time.dat'
       OPEN(7, FILE = CTemp, STATUS = 'NEW')
       CLOSE(7)

       CTemp = 'TempAve.dat'
       OPEN(7, FILE = CTemp, STATUS = 'NEW')
       CLOSE(7)
    ELSE IF( SFLAG .EQ. 1 ) THEN
        t = t+dt
    END IF

    PRINT *, '==============================='
    PRINT *, 'START MAIN PROGRAM'

    mainloop: DO it = Startit, Maxit

       CALL EvolveData(&
            &Nr, Nth, Nphi, Mr, Lmax, Lgrid,&
            &GLQWeights, GLQZeros,&
            &rmin, rmax,&
            &r, rho, theta, phi,&
            &epsStar, alpStar,&
            &SolPhi,&
            &t, dt,&
            &a)

       CALL GetResults(&
            &Nr, Nth, Nphi, Mr, Lmax, Lgrid, SpM,&
            &GLQWeights, GLQZeros,&
            &r, rho, theta, phi,&
            &a,&
            &it, WriteSit,&
            &thetaSp, phiSp,&
            &TempAve) 

       !--------------------------------------------------------!
       !     Writing OUTPUTS into files                         !
       !--------------------------------------------------------!

       IF( MOD(it,Writeait) .EQ. 0 .OR. (it .EQ. Maxit) ) &
         & CALL Writea(Mr+1, 2, Lmax+1, Lmax+1, a, it)
       
       CTemp = 'Time.dat'
       OPEN(7, FILE = CTemp, ACCESS = 'APPEND', STATUS = 'OLD')
       WRITE(7,*) t
       CLOSE(7)
 
       CTemp = 'TempAve.dat'
       OPEN(7, FILE = CTemp, ACCESS = 'APPEND', STATUS = 'OLD')
       WRITE(7,*) TempAve
       CLOSE(7)

       PRINT *, 'Iteration #:', it
       PRINT *, 'Time= ', t
       PRINT *, 'Average Temperature= ', TempAve
       PRINT *, '------------------'

       IF( ((t .LT. tfinal) .AND. (tdir .LT. 0.0D0)) .OR.&
          &((t .GT. tfinal) .AND. (tdir .GT. 0.0D0)) ) THEN
           CALL Writea(Mr+1, 2, Lmax+1, Lmax+1, a, it)
           EXIT
       END IF

    END DO mainloop

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
  END PROGRAM SHF

