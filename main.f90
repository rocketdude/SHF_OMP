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

    INTEGER*4, PARAMETER ::        Mr       = 40
    INTEGER*4, PARAMETER ::        Lmax     = 64
    INTEGER*4, PARAMETER ::        Lgrid    = 64
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

    INTEGER*4           WriteUit    
    !iteration at which we output U into file
    INTEGER*4           Writeg_rrUsqrdit 
    !Iteration at which we output gRR*U^2
    INTEGER*4           WriteSit        
    !iteration at which we output S into file
    INTEGER*4           Writeait        
    !iteration at which we output the coef. a_nlm

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
    REAL*8            rootsign    !Choose the sign of the root (+/-)

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

    REAL*8              U(Nth,Nphi) !Radial distance of the light cone
    REAL*8              Uave        !Average radial distance of the light cone
    REAL*8              g_rrUsqrd(Nth,Nphi)!g_rr * U^2
    REAL*8              g_rrUsqrdAve  !Average g_rr * U^2
    REAL*8              USp(SpM)      !The values of U in specified directions

    INTEGER*4           n           !Degree of Chebyshev polynomial
    INTEGER*4           l, ml       !Degree of spherical harmonics
    INTEGER*4           reinit      !Iteration at which we reinitialize

!--------------------------------------------------------!
!     Declare Metric Data                                !
!--------------------------------------------------------!

    ! Values of the upper index metric

    ! Spatial metric (gamma's)--all upper indices
    ! Remember that gamma^ij is the inverse of the SPATIAL metric alone
    REAL*8              gRR(Nr,Nth,Nphi)
    REAL*8              gThTh(Nr,Nth,Nphi)
    REAL*8              gPhiPhi(Nr,Nth,Nphi)
    REAL*8              gRTh(Nr,Nth,Nphi)      
    REAL*8              gRPhi(Nr,Nth,Nphi)
    REAL*8              gThPhi(Nr,Nth,Nphi)

    REAL*8              alpha(Nr,Nth,Nphi)    !Lapse function
    REAL*8              betaR(Nr,Nth,Nphi)    !Shift function beta^{r}
    REAL*8              betaTh(Nr,Nth,Nphi)   !Shift function beta^{theta}
    REAL*8              betaPhi(Nr,Nth,Nphi)  !Shift function beta^{phi}

    ! Values of the upper index metric at different times
    REAL*8              BgRR(TP,Nr,Nth,Nphi)
    REAL*8              BgThTh(TP,Nr,Nth,Nphi)
    REAL*8              BgPhiPhi(TP,Nr,Nth,Nphi)
    REAL*8              BgRTh(TP,Nr,Nth,Nphi)
    REAL*8              BgRPhi(TP,Nr,Nth,Nphi)
    REAL*8              BgThPhi(TP,Nr,Nth,Nphi)

    REAL*8              Balpha(TP,Nr,Nth,Nphi)
    REAL*8              BbetaR(TP,Nr,Nth,Nphi)
    REAL*8              BbetaTh(TP,Nr,Nth,Nphi)
    REAL*8              BbetaPhi(TP,Nr,Nth,Nphi)

    ! Schwarzschild metric parameters
    REAL*8            Mass   !Mass of Schwarzschild black hole (for testing)

    ! Parameters to read metric from HDF5 files
    INTEGER(HSIZE_T)  bufsize(3)    !Buffer size to read the metric from HDF5
    INTEGER*4         nchunks       !Number of chunks that output the metric
    INTEGER*4         it_data(TP)   !The Einstein Toolkit iteration value
    INTEGER*4         it_data_max   !Maximum Einstein Toolkit iteration value
    INTEGER*4         it_data_min   !Minimum Einstein Toolkit iteration value
    INTEGER*4         delta_it_data !Difference in it_data between metric values
    INTEGER*4         readdata      !Index of it_data that needs to be read
    INTEGER*4         it_data_test

    REAL*8            t_data(TP)    !Data times
    REAL*8            t_thresh      !Threshold time, beyond this read new data
    INTEGER*4         Maxit_allowed !The amount of iterations that we can have
    
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
    t = 50.0D0                      !Last time from previous run
    Startit = 6751                  !Startit = last iteration + 1
    aFile = 'a10.dat'

    !Termination conditions
    Maxit = 20000
    tfinal = 5.0D0

    IF( SFLAG .EQ. 0 ) THEN
       Startit = 1                  !Starting from iteration 1
    END IF

    X0 = 0.700D0                    !Initial axes of the spheroids in the
    Y0 = 0.700D0                    !X, Y, and Z directions
    Z0 = 0.600D0

    !Simulation parameters                              
    !Note: negative rootsign, positive lapse & shift functions, 
    !      and negative tdir give EH finder
    eps =  2.22044604925031308D-016 !Machine epsilon (precalculate)
    c = 0.1D0
    cfl = 1.508D0                   !Depends on which SSP-Runge-Kutta used
                                    !SSPRK(5,4)=>cfl=1.508 and 
                                    !SSPRK(3,3)=>cfl=1.0
    rootsign = -1.0D0               !Choose the root sign, 
                                    !either 1.0D0 or -1.0D0 (depends on metric)
                                    !e.g. for alpha=(+), beta=(+), dt=(-), 
                                    !     -1.0D0 is an event horizon finder
    tdir = -1.0D0                   !Direction of time, choose +1.0D0 or -1.0D0
    reinit = 15

    !Spherical grid parameters
    rmax = 1.20D0                   !maximum value of r
    rmin = 0.20D0                   !minimum value of r

    !Additional directions we'd like to compute U
    thetaSp = (/ 0.0D0, PI, PI/2.0D0, PI/2.0D0, PI/2.0D0, PI/2.0D0 /)
    phiSp = (/ 0.0D0, 0.0D0, 0.0D0, PI/2.0D0, PI, 1.5D0*PI /)
    
    !Parameters related to reading HDF5 files--do h5dump to check these
    nchunks = 16
    bufsize(1) = 50 !Buffer sizes need to be bigger than datasets
    bufsize(2) = 50
    bufsize(3) = 50
    it_data_max = 8000
    it_data_min = 0
    delta_it_data = 4

    !Schwarzschild metric parameter (preliminary tests only)
!!$    Mass = 0.45D0                   !Only used with Schwarzschild metric

!--------------------------------------------------------!
!     Output Parameters                                  !
!--------------------------------------------------------!

    WriteUit         = 1000
    Writeg_rrUsqrdit = 1000
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

    !Set the number of MKL threads (dynamic is automatic)
    CALL MKL_SET_DYNAMIC(.TRUE.)

    !Inquire number of threads
    nthreads = omp_get_max_threads()
    WRITE(*,*) 'No. of threads = ', nthreads

    !Set up Chebyshev roots and DH collocation points
    !$OMP PARALLEL DO
    DO i = 0, Mr
       rho(i+1) = -COS(PI*DBLE(2*i + 1) / DBLE(2*(Mr+1)) )
       r(i+1) = 0.5D0*( (rmax + rmin) + (rmax - rmin)*rho(i+1) )
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
    PRINT *, 'Number of metric data chunks =', nchunks

!--------------------------------------------------------!
!     Find smallest spatial increment                    !
!--------------------------------------------------------!

    dx = rmax - rmin
    DO i = 1, Mr
       IF( (r(i+1) - r(i)) .LT. dx ) THEN
          dx = r(i+1) - r(i)
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
            & r, rho, theta, phi,&
            & a)
    END IF

!--------------------------------------------------------!
!     Evaluate Metric                                    !
!--------------------------------------------------------!
    !This Schwarzschild metric is only for preliminary calculations
    !purposes only
!!$    CALL GetSchwarzschildMetric(&
!!$         &Mass,&
!!$         &M, Mr, NP,&
!!$         &r, theta, phi,&
!!$         &alpha,&
!!$         &betaR, betaTh, betaPhi,&
!!$         &gRR, gThTh, gPhiPhi,&
!!$         &gRTh, gRPhi, gThPhi)

    PRINT *, '==============================='


!--------------------------------------------------------!
!     Evolve Eikonal S                                   !
!--------------------------------------------------------!

    IF( SFLAG .EQ. 0 ) THEN
       !Initialize output files and write values for iteration 0
       CTemp = 'Time.dat'
       OPEN(7, FILE = CTemp, STATUS = 'NEW')
       CLOSE(7)
       CTemp = 'Uave.dat'
       OPEN(7, FILE = CTemp, STATUS = 'NEW')
       CLOSE(7)
       CTemp = 'USp.dat' 
       OPEN(7, FILE = CTemp, STATUS = 'NEW')
       CLOSE(7)
       CTemp = 'g_rrUsqrdAve.dat'
       OPEN(7, FILE = CTemp, STATUS = 'NEW')
       CLOSE(7)
    ELSE IF( SFLAG .EQ. 1 ) THEN
        t = t+dt
    END IF

    !Initializes read metric logics
    DO i = 1, TP
        it_data(i) = it_data_max - (i-1)*delta_it_data
    END DO

    readdata = 0
    t_thresh = tdir * ABS(1.0D30)

    PRINT *, '==============================='
    PRINT *, 'START MAIN PROGRAM'

    mainloop: DO it = Startit, Maxit
       
       !--------------------------------------------------------!
       !     Read Metric and Evolve Data                        !
       !--------------------------------------------------------!

        IF ( ((tdir .LT. 0.0D0) .AND. (t .LT. t_thresh)) .OR. &
            &((tdir .GT. 0.0D0) .AND. (t .GT. t_thresh)) ) THEN
            readdata = MAXLOC(it_data,1)
            it_data_test = MINVAL(it_data) - delta_it_data

            IF (it_data_test .LT. it_data_min) THEN
                readdata = -1
            ELSE
                it_data(readdata) = it_data_test
            END IF
        END IF
                
        CALL GetMetricAtCurrentTime(&
             &Nr, Nth, Nphi, TP,&
             &readdata,&
             &t, t_thresh, tdir,&
             &t_data, it_data,&    
             &nchunks,&
             &bufsize,&
             &r, theta, phi,&
             &Balpha, BbetaR, BbetaTh, BbetaPhi,&
             &BgRR, BgThTh, BgPhiPhi,&
             &BgRTh, BgRPhi, BgThPhi,&
             &alpha, betaR, betaTh, betaPhi,&
             &gRR, gThTh, gPhiPhi,&
             &gRTh, gRPhi, gThPhi)

       CALL EvolveData(&
            &Nr, Nth, Nphi, Mr, Lmax, Lgrid,&
            &GLQWeights, GLQZeros,&
            &rootsign, rmin, rmax,&
            &rho, theta, phi,&
            &alpha,&
            &betaR, betaTh, betaPhi,&
            &gRR, gThTh, gPhiPhi,&
            &gRTh, gRPhi, gThPhi,&
            &t, dt,&
            &a)
 
       CALL FindU(&
            &Nr, Nth, Nphi, Mr, Lmax, Lgrid, SpM,&
            &GLQWeights, GLQZeros,&
            &gRR, gThTh, gPhiPhi,&
            &gRTh, gRPhi, gThPhi,&
            &r, rho, theta, phi,&
            &a,&
            &it, WriteSit,&
            &U, Uave, USp,&
            &thetaSp, phiSp,&
            &g_rrUsqrd, g_rrUsqrdAve)


       !--------------------------------------------------------!
       !     Writing OUTPUTS into files                         !
       !--------------------------------------------------------!

       IF( MOD(it,WriteUit) .EQ. 0 ) CALL WriteU(Nth, Nphi, U, it)
       IF( MOD(it,Writeg_rrUsqrdit) .EQ. 0 ) &
         & CALL WriteGRRUU(Nth, Nphi, g_rrUsqrd, it)
       IF( MOD(it,Writeait) .EQ. 0 .OR. (it .EQ. Maxit) ) &
         & CALL Writea(Mr+1, 2, Lmax+1, Lmax+1, a, it)
       
       CTemp = 'Time.dat'
       OPEN(7, FILE = CTemp, ACCESS = 'APPEND', STATUS = 'OLD')
       WRITE(7,*) t
       CLOSE(7)

       CTemp = 'Uave.dat'
       OPEN(7, FILE = CTemp, ACCESS = 'APPEND', STATUS = 'OLD')
       WRITE(7,*) Uave
       CLOSE(7)

       CTemp = 'USp.dat'
       OPEN(7, FILE = CTemp, ACCESS = 'APPEND', STATUS = 'OLD')
       DO i = 1,SpM
            WRITE(7,*) USp(i)
       END DO
       CLOSE(7)
    
       CTemp = 'g_rrUsqrdAve.dat'
       OPEN(7, FILE = CTemp, ACCESS = 'APPEND', STATUS = 'OLD')
       WRITE(7,*) g_rrUsqrdAve
       CLOSE(7)
    
       !--------------------------------------------------------!
       !     SMOOTHING: Reinitializing S                        !
       !--------------------------------------------------------!

       IF( MOD(it, reinit) .EQ. 0 ) THEN
           
           CALL ReinitializeData(&
                   &Nr, Nth, Nphi,&
                   &Mr, Lmax, Lgrid,&
                   &GLQWeights, GLQZeros,&
                   &r, rho, theta, phi,&
                   &c, U, a)

       END IF

       PRINT *, 'Iteration #:', it
       PRINT *, 'Time= ', t
       PRINT *, 'Average R= ', Uave
       PRINT *, 'Metric data used: '
       DO i = 1, TP
           PRINT *, 'Iteration#:', it_data(i)
       END DO
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

