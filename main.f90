!======================================================! 
! Spectral Horizon Finder (SHF)                        ! 
!======================================================!

! This program is an event horizon finder using the pseudospectral methods.
! The eikonal data S(t, r, theta, phi) is represented as orthogonal functions 
! of Chebyshev polynomials & spherical harmonics.
! In other words, we are approximating S(t, r, theta, phi) ~= a_nlm(t) T_n (rho) Y_lm (theta, phi).


  PROGRAM              SHF
    USE                omp_lib
    USE                HDF5
    !USE                H5LT

    IMPLICIT           none

    INTEGER*4, PARAMETER ::        Mr = 100
    INTEGER*4, PARAMETER ::        M = 8
    INTEGER*4, PARAMETER ::        NP = (Mr+1)*(M*M)

    ! Mr is the degree of the radial Chebyshev polynomial
    ! M is the degree of angular spherical harmonics with maximum L = M-1
    ! M has to be multiples of two
    ! NP is the number of coefficients, 4*NP is the number of collocation points

    REAL*8, PARAMETER ::         PI = 3.141592653589793238462643383279502884197D0

!-------------------------------------------------------!
!     Declare Commons                                   !
!-------------------------------------------------------!

!-------------------------------------------------------!
!     DATA I/O                                          !
!-------------------------------------------------------!

    CHARACTER*32      CTemp
    CHARACTER*32      aFile
    LOGICAL           FileExist

    INTEGER*4           WriteUit !iteration at which we output U into file
    INTEGER*4           WriteSit !iteration at which we output S into file
    INTEGER*4           Writeait !iteration at which we output the coefficients a_nlm into file

!--------------------------------------------------------!
!     Declare Numerical Inputs                           ! 
!--------------------------------------------------------!

    INTEGER*4           it     !Iteration counter
    INTEGER*4           Startit!Starting iteration
    INTEGER*4           Maxit  !Maximum iteration

    INTEGER*4           SFLAG  !If SFLAG = 1, then we are continuing previous run
    INTEGER*4           LWORK  !Size of the WORK matrix, used for inverting matrices
    INTEGER*4           IFLAG  !IFLAG is inverse flag. IFLAG = 1, AFinv.dat is available

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
    ! Here we note that the coordinate locations do not correspond to the collocation points
    ! We only use half of the theta coordinates for the collocation points

    REAL*8            t           !Current time
    REAL*8            tdir        !Direction of time (+1 for forward, -1 for backward)
    REAL*8            dt          !Time increment
    REAL*8            dx          !Smallest distance between two points (used to determine dt)
    REAL*8            rdtheta     !Smallest distance in the theta direction
    REAL*8            rsinthdphi  !Smallest distance in the phi direction
    REAL*8            cputime_start, cputime_stop    !used for calculating the program speed
    REAL*8            eps         !Machine epsilon (needs to be calculated before running this)
    REAL*8            rootsign    !Choose the sign of the root (+/-)

    REAL*8            rho(Mr+1)   !Canonical radial coordinate for Chebyshev polynomials
    REAL*8            r(Mr+1)     !The radial coordinate r according to the usual sph. coord. system
    REAL*8            rmax, rmin  !rmax and rmin define the radial domain in this program,
                                  !used to calculate the canonical radial coordinate rho

    REAL*8            theta(2*M)  !The angle theta according to the usual sph. coordinate system
    REAL*8            phi(2*M)    !The angle phi according to the usual sph. coordinate system

    REAL*8            R0          !Initial radius of the light 'cone'

    REAL*8            U(2*M,2*M)  !Radial distance of the light cone in each theta and phi direction
    REAL*8            Uave        !Average radial distance of the light cone at a specific time step

    COMPLEX*16        a(NP)       !Values of the time-dependent coefficients a_nlm(t)
                                  !vectorized in loop n:l:m
    COMPLEX*16        AF(4*NP, NP) !Matrix of Tn(r)Ylm(theta, phi)
    COMPLEX*16        AFinv(NP, 4*NP) !Inverse of matrix AF
    COMPLEX*16        B(NP, 4*NP) !Matrix of quadrature
    COMPLEX*16        Dth(4*NP, NP) !Matrix used to take derivative of S w.r.t. theta
    COMPLEX*16        Dphi(4*NP, NP)!Matrix used to take derivative of S w.r.t phi
    REAL*8            F(NP)       !Exponential filter for the Chebyshev polynomials

    INTEGER*4           Lmax        !Maximum degree of l of the spherical harmonics
    INTEGER*4           n           !Degree of Chebyshev polynomial
    INTEGER*4           l, ml       !Degree of spherical harmonics
    INTEGER*4           reinit      !Iteration at which we smooth the data by reinitializing

!--------------------------------------------------------!
!     Declare Metric Data                                !
!--------------------------------------------------------!

    ! Values of the upper index metric at the different collocation points
    ! The collocation points are (r, theta, phi)

    ! Spatial metric (gamma's)--all upper indices
    ! Remember that gamma^ij is the inverse of the SPATIAL metric alone
    REAL*8        gRR(4*NP)
    REAL*8        gThTh(4*NP)
    REAL*8        gPhiPhi(4*NP)
    REAL*8        gRTh(4*NP)      
    REAL*8        gRPhi(4*NP)
    REAL*8        gThPhi(4*NP)

    REAL*8        alpha(4*NP)    !Lapse function
    REAL*8        betaR(4*NP)    !Shift function beta^{r}
    REAL*8        betaTh(4*NP)   !Shift function beta^{theta}
    REAL*8        betaPhi(4*NP)  !Shift function beta^{phi}

    ! Schwarzschild metric parameters
    REAL*8            Mass   !Mass of Schwarzschild black hole located at the center

    ! Parameters to read metric from HDF5 files
    INTEGER(HSIZE_T)  bufsize(3)  !Buffer size to read the metric from HDF5 files
    INTEGER*4         iterratio   !Ratio of the iteration values in the EinsteinToolkit compared to this code's iteration
    INTEGER*4         nchunks     !Number of chunks that the metric is divided into
    INTEGER*4         it_data     !The iteration value for the metric data
    REAL*8            dt_data     !The value of dt for the metric data

!--------------------------------------------------------!
!     Declare Local Parameters                           !
!--------------------------------------------------------!

    INTEGER*4           i,j,k       !Counters for DO loops
    INTEGER*4           crow        !Row counter
    INTEGER*4           nthreads    !Number of threads

!========================================================!
!     MAIN PROGRAM                                       !
!========================================================!

    Lmax = M - 1
!--------------------------------------------------------!
!     Parameters                                         !
!--------------------------------------------------------!
    
    !Parameters related to initial conditions
    SFLAG = 0                       !If SFLAG = 1, we are continuing previous run: change t, Startit and aFile
    t = -8.1617491479814266D0       !Last time from previous run
    Startit = 6751                  !Startit = last iteration + 1
    Maxit = 500
    aFile = 'a10.dat'

    IF( SFLAG .EQ. 0 ) THEN
       t  = 0.0D0                   !Initial time
       Startit = 1                  !Starting from iteration 1
    END IF
    R0 = 0.52D0                     !Initial radius of the null surface

    !Simulation parameters                              
    !Note: negative rootsign, positive lapse & shift functions, and negative tdir give EH finder
    eps =  2.22044604925031308D-016 !Machine epsilon (calculate everytime you change machine)
    c = 0.1D0
    cfl = 1.508D0
    dt = 2.50D-3                    !Just enter the absolute value, change the sign by changing tdir
    rootsign = -1.0D0               !Choose the root sign, either 1.0D0 or -1.0D0 (depends on the metric)
    tdir = -1.0D0                   !Direction of time, choose +1.0D0 or -1.0D0
    reinit = 15

    IFLAG = 1                       !IFLAG = 0 ==> calculate AFinv (slow)
                                    !IFLAG = 1 ==> read AFinv from AFinv.dat
                                    !IFLAG = -1 ==> calculate AFinv and write it into AFinv.dat
    !Note:LWORK needs to be changed everytime we change M or Mr, AFinv has to be recalculated
    LWORK = 42617152                !Optimized size of the WORK matrix, Put LWORK = -1 to query for the optimal size

    !Spherical grid parameters
    rmax = 1.10D0                   !maximum value of r
    rmin = 0.0D0                    !minimum value of r

    !Parameters related to reading HDF5 files--do h5dump to check these
    iterratio = 4
    nchunks = 12
    bufsize(1) = 50                 !Buffer size; has to be bigger than the size of each dataset
    bufsize(2) = 50
    bufsize(3) = 50
    dt_data = 5.00D-3               !Size of time step specified by the data
    !dt = 0.0125 but rl=8 is only updated once every four iterations, so dt = 0.050

    !Schwarzschild metric parameter (preliminary tests only)
!!$    Mass = 0.45D0                   !Only used with Schwarzschild metric

!--------------------------------------------------------!
!     Output Parameters                                  !
!--------------------------------------------------------!

    WriteUit = 10
    WriteSit = 100
    Writeait = 100

!--------------------------------------------------------!
!     Echo certain parameters                            !
!--------------------------------------------------------!

    IF( (MOD(M,2) .EQ. 0) .OR. (M .EQ. 1) ) THEN
       PRINT *, 'M =', M
       PRINT *, 'Maximum L =', Lmax
    ELSE
       PRINT *, 'M has to be a multiple of two (or M=1)'
       STOP
    END IF

!!$    PRINT *, 'Mass of blackhole =', Mass
    PRINT *, '# of iterations =', (Maxit-Startit+1)
    PRINT *, 'Reinitializing every ', reinit, 'iterations'

!--------------------------------------------------------!
!     Timer Start                                        !
!--------------------------------------------------------!	  

    CALL cpu_time( cputime_start )

!--------------------------------------------------------!
!     Creating mesh & Inquiring Threads                  !
!--------------------------------------------------------!

    !Set the number of MKL threads (dynamic is automatic)
    CALL MKL_SET_DYNAMIC(.TRUE.)

    !Inquire number of threads
    nthreads = omp_get_max_threads()
    WRITE(*,*) 'No. of threads = ', nthreads

    !$OMP PARALLEL
    !$OMP DO
    DO i = 0, Mr
       rho(i+1) = -COS(PI*DBLE(2*i + 1) / DBLE(2*(Mr+1)) )
       r(i+1) = 0.5D0*( (rmax + rmin) + (rmax - rmin)*rho(i+1) )
    END DO
    !$OMP END DO

    !$OMP DO
    DO j = 0, (2*M-1)
       theta(j+1) = PI*(DBLE(j) + 0.5D0 )/DBLE(2*M)
       phi(j+1) = PI*2.0D0*(DBLE(j) + 0.5D0)/DBLE(2*M)
    END DO
    !$OMP END DO
    !$OMP END PARALLEL

!--------------------------------------------------------!
!     Find smallest spatial increment                    !
!--------------------------------------------------------!

    dx = rmax - rmin
    DO i = 1, Mr
       IF( (r(i+1) - r(i)) .LT. dx ) THEN
          dx = r(i+1) - r(i)
       END IF
    END DO

!!$    rdtheta = MINVAL(r)*PI/DBLE(2*M)
!!$    rsinthdphi = MINVAL(r)*MINVAL(SIN(theta))*PI/DBLE(M)
!!$
!!$    IF( rdtheta .LT. dx ) THEN
!!$       dx = rdtheta
!!$    ELSEIF( rsinthdphi .LT. dx ) THEN
!!$       dx = rsinthdphi
!!$    END IF

    IF( ABS(dt) .GT. ABS(tdir * cfl * dx) ) THEN
        PRINT *, 'Warning: dt might be greater than Courant limit!'
    END IF
    dt = ABS(dt) * tdir

!--------------------------------------------------------!
!     Initial Data                                       !
!--------------------------------------------------------!

    PRINT *, '==============================='
    PRINT *, 'INITIALIZING PROGRAM'

    CALL GetMatrices(& 
         & M, Mr, NP, Lmax, LWORK, IFLAG,&
         & r, rho, theta, phi,&
         & AF, AFinv, B, Dth, Dphi)

    IF( SFLAG .EQ. 1 ) THEN
       PRINT *, 'Continuing run'
       PRINT *, 'Previous iteration=' ,(Startit-1)

       CALL Read1dC(NP, a, aFile)

    ELSE

       CALL GetInitialData(& 
            & M, Mr, NP,&
            & c,&
            & R0,&
            & r, theta, phi,&
            & B,&
            & a)

    END IF

    CALL GetFilter(& 
         & M, Mr, NP, Lmax,&
         & eps,&
         & F)

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

    IF( Startit .EQ. 1 ) THEN
       !Initialize output files and write values for iteration 0
       CTemp = 'Time.dat'
       OPEN(7, FILE = CTemp, STATUS = 'NEW')
       WRITE(7,*) t
       CLOSE(7)
       CTemp = 'Uave.dat'
       OPEN(7, FILE = CTemp, STATUS = 'NEW')
       WRITE(7,*) R0
       CLOSE(7)
    END IF

    PRINT *, '==============================='
    PRINT *, 'START MAIN LOOP'

    mainloop: DO it = Startit, Maxit
       
       t = t+dt
       
       !--------------------------------------------------------!
       !     Read Metric and Evolve Data                        !
       !--------------------------------------------------------!

       it_data = iterratio * (Maxit - it)

       CALL GetMetric(&
            &M, Mr, NP,&
            &dt, dt_data,&
            &it_data, nchunks,&
            &bufsize,&
            &r, theta, phi,&
            &alpha,&
            &betaR, betaTh, betaPhi,&
            &gRR, gThTh, gPhiPhi,&
            &gRTh, gRPhi, gThPhi)


       CALL EvolveData(&
            &M, Mr, NP, Lmax,&
            &rootsign, rmin, rmax,&
            &alpha,&
            &betaR, betaTh, betaPhi,&
            &gRR, gThTh, gPhiPhi,&
            &gRTh, gRPhi, gThPhi,&
            &dt,&
            &AF, AFinv, Dth, Dphi, F,&
            &a)

       CALL FindU(&
            &M, Mr, NP,&
            &r,&
            &AF,a,&
            &it, WriteSit,&
            &U, Uave)


       !--------------------------------------------------------!
       !     Writing OUTPUTS into files                         !
       !--------------------------------------------------------!

       IF( MOD(it,WriteUit) .EQ. 0 ) CALL WriteU(2*M, 2*M, U, it)
       IF( MOD(it,Writeait) .EQ. 0 .OR. (it .EQ. Maxit) ) CALL Writea(NP, a, it)
       
       CTemp = 'Time.dat'
       OPEN(7, FILE = CTemp, ACCESS = 'APPEND', STATUS = 'OLD')
       WRITE(7,*) t
       CLOSE(7)

       CTemp = 'Uave.dat'
       OPEN(7, FILE = CTemp, ACCESS = 'APPEND', STATUS = 'OLD')
       WRITE(7,*) Uave
       CLOSE(7)

       !--------------------------------------------------------!
       !     SMOOTHING: Reinitializing S                        !
       !--------------------------------------------------------!

       IF( MOD(it, reinit) .EQ. 0 ) THEN

          CALL ReinitializeData(& 
               & M, Mr, NP,&
               & c,&
               & U,&
               & r, theta, phi,&
               & B,&
               & a)

       END IF

       PRINT *, 'Iteration #:', it
       PRINT *, 'Time= ', t
       PRINT *, 'Average R= ', Uave
       PRINT *, '------------------'

    END DO mainloop

    PRINT *, 'END PROGRAM'
    PRINT *, '==============================='

!--------------------------------------------------------!
!     Timer Stop                                         !
!--------------------------------------------------------!

    CALL cpu_time( cputime_stop )
    WRITE(*,*) 'Elapsed CPU Time = ', cputime_stop - cputime_start

    STOP
  END PROGRAM SHF

