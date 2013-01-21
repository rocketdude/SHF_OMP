!---------------------------------------------------------!
!   Solve the equation using globally convergent          !
!   Newton-Raphson method                                 !
!---------------------------------------------------------!

    SUBROUTINE  SolveEquation(&
                &Nth, Nphi, Lmax, Lgrid, GLQWeights, GLQZeros,&
                &R, theta, phi,&
                &Maxit, tolF, tolMin,&
                &LWORK,&
                &a)

        USE         omp_lib
        IMPLICIT    none
    
        !Declare passed variables
        INTEGER*4               Nth, Nphi, Lmax, Lgrid, Maxit
        INTEGER*4               LWORK
        REAL*8                  GLQWeights(Lgrid+1), GLQZeros(Lgrid+1)
        REAL*8                  R, theta(Nth), phi(Nphi)
        REAL*8                  tolF,tolMin
        REAL*8                  a(2,Lmax+1,Lmax+1)

        !Declare local variables
        INTEGER*4               i, j, k, it
        REAL*8                  errF, errA

        REAL*8                  aV((Lmax+1)**2),aVOld((Lmax+1)**2)
        REAL*8                  Fn((Lmax+1)**2)
        REAL*8                  g((Lmax+1)**2), p((Lmax+1)**2)
        REAL*8                  Jacobian((Lmax+1)**2,(Lmax+1)**2)
        REAL*8                  invJacobian((Lmax+1)**2,(Lmax+1)**2)

        REAL*8                  fOld,f,stpmax
        REAL*8                  tolX,stpmx
        LOGICAL                 check

        !Needed for solving the J . dx = -F
        REAL*8                  WORK(MAX(1,LWORK))
        INTEGER*4               RANK,INFO
        REAL*8                  S((Lmax+1)**2)

        !Main subroutine

        CALL SHCilmToVector(a,aV,Lmax)
        tolX = epsilon(aV)
        stpmx = 100.0D0

        CALL RadiationFunction(Nth,Nphi,Lmax,Lgrid,&
            &GLQWeights,GLQZeros,R,theta,phi,&
            &aV,Fn)

        f = 0.5D0*DOT_PRODUCT(Fn,Fn)

        IF( MAXVAL(ABS(Fn)) < 0.01D0*tolF ) THEN
            check = .false. !check if initial guess is a solution
            PRINT *, 'Initial guess is a solution'          
            RETURN
        END IF

        stpmax = stpmx*MAX(SQRT(DOT_PRODUCT(aV,aV)),DBLE((Lmax+1)**2))

        DO it=1,Maxit
        
            !Compute the Jacobian by numerical differentiation
            CALL FDJacobian(Nth,Nphi,Lmax,Lgrid,&
                &GLQWeights,GLQZeros,R,theta,phi,&
                &aV,Jacobian)

            !Check the maximum value of the Jacobian vs. the minimum value
            PRINT *, 'Maximum value = ', MAXVAL(Jacobian)
            PRINT *, 'Minimum value = ', MINVAL(Jacobian)
            PRINT *, 'Ratio of max/min = ', MAXVAL(Jacobian)/MINVAL(Jacobian)

            g = MATMUL(Fn(:),Jacobian(:,:))
            aVOld = aV
            fold = f
            p = -1.0D0*Fn

            !Compute da by using SVD
            CALL DGELSS((Lmax+1)**2,(Lmax+1)**2,1,Jacobian,&
                    &(Lmax+1)**2,p,(Lmax+1)**2,S,-1.0D0,RANK,WORK,LWORK,INFO)

            IF( LWORK .EQ. -1 ) THEN
                PRINT *, 'Work space query, optimal LWORK =', WORK(1)
                STOP
            END IF

            IF( INFO .GT. 0 ) THEN
                PRINT *, 'INFO = ', INFO, ' Failed to converge'
                STOP
            ELSEIF( INFO .LT. 0 ) THEN
                PRINT *, 'INFO = ', INFO, ' Illegal value'
                STOP
            END IF 

            !Do linesearch
            CALL lnsrch(aVOld,fOld,g,p,aV,f,stpmax,check,Fn,&
                      &Nth,Nphi,Lmax,Lgrid,GLQWeights,GLQZeros,R,theta,phi)

            !Check if residual is small enough
            errF = MAXVAL( ABS(Fn) )
            !errF = SQRT(DOT_PRODUCT(Fn,Fn)) / (Lmax+1)**2
            IF( errF .LE. tolF ) THEN
                PRINT *, 'Solved: error in residual is within tolerance'
                CALL SHVectorToCilm(aV,a,Lmax)
                RETURN
            END IF

!!$            IF( check ) THEN !Check for gradient of f to be 0
!!$                check=(MAXVAL(ABS(g)*MAX(ABS(aV),1.0D0) / &
!!$                        &MAX(Fn,0.5D0*(Lmax+1)**2)) < tolMin)
!!$                PRINT *, 'Gradient is 0 but solution had yet to converge'
!!$                CALL SHVectorToCilm(aV,a,Lmax)
!!$                RETURN
!!$            END IF

            errA = MAXVAL(ABS(aV-aVOld)/MAX(ABS(aV),1.0D0))
            IF( errA < tolX ) THEN
                PRINT *, 'delta a converges'
                CALL SHVectorToCilm(aV,a,Lmax)
                RETURN
            END IF

            !Output to user
            PRINT *, 'Iteration # =', it
            PRINT *, 'error in residual =', errF
            PRINT *, 'error in solution =', errA
            PRINT *, '========================='
            
        END DO

        PRINT *, 'ERROR: Maximum number of iterations reached'
        CALL SHVectorToCilm(aV,a,Lmax)
        RETURN
    END SUBROUTINE


!=================================================================!

    SUBROUTINE lnsrch(xold,fold,g,p,x,f,stpmax,check,Fn,&
                      &Nth,Nphi,Lmax,Lgrid,GLQWeights,GLQZeros,R,theta,phi)

    IMPLICIT NONE

    !Needed to evaluate the radiation function
    INTEGER*4                   :: Nth,Nphi,Lmax,Lgrid
    REAL*8                      :: GLQWeights(Lgrid+1),GLQZeros(Lgrid+1)
    REAL*8                      :: R,theta(Nth),phi(Nphi)

    REAL*8, INTENT(IN)          :: xold((Lmax+1)**2),g((Lmax+1)**2)
    REAL*8, INTENT(IN)          :: fold,stpmax
    REAL*8, INTENT(INOUT)       :: p((Lmax+1)**2)

    REAL*8, INTENT(OUT)         :: x((Lmax+1)**2)
    REAL*8, INTENT(OUT)         :: f
    LOGICAL, INTENT(OUT)        :: check
    REAL*8, INTENT(OUT)         :: Fn((Lmax+1)**2)

    !Local parameters
    REAL*8, PARAMETER           :: ALF=1.0D-4, TOLX=epsilon(x)
    INTEGER*4                   :: ndum
    REAL*8 :: a,alam,alam2,alamin,b,disc,f2,pabs,rhs1,rhs2,slope,tmplam

    !Main subroutine
    check=.false.
    pabs=sqrt(dot_product(p,p))
    if (pabs > stpmax) p(:)=p(:)*stpmax/pabs
    slope=dot_product(g,p)
    if (slope >= 0.0d0) STOP "***roundoff problem in lnsrch***"

    alamin=TOLX/maxval(abs(p(:))/max(abs(xold(:)),1.0d0))
    alam=1.0d0
    do
        x(:)=xold(:)+alam*p(:)
       
        !We are trying to minimize f = 1/2*F (dot) F 
        CALL RadiationFunction(Nth,Nphi,Lmax,Lgrid,&
            &GLQWeights,GLQZeros,R,theta,phi,&
            &x,Fn)
        f = 0.5D0*dot_product(Fn,Fn)

        if (alam < alamin) then
            x(:)=xold(:)
            check=.true.
            RETURN
        else if (f <= fold+ALF*alam*slope) then
            RETURN
        else
            if (alam == 1.0d0) then
                tmplam=-slope/(2.0d0*(f-fold-slope))
            else
                rhs1=f-fold-alam*slope
                rhs2=f2-fold-alam2*slope
                a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
                b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/&
                    (alam-alam2)
                if (a == 0.0d0) then
                    tmplam=-slope/(2.0d0*b)
                else
                    disc=b*b-3.0d0*a*slope
                    if (disc < 0.0d0) then
                        tmplam=0.5d0*alam
                    else if (b <= 0.0d0) then
                        tmplam=(-b+sqrt(disc))/(3.0d0*a)
                    else
                        tmplam=-slope/(b+sqrt(disc))
                    end if
                end if
                if (tmplam > 0.5d0*alam) tmplam=0.5d0*alam
            end if
        end if
        alam2=alam
        f2=f
        alam=max(tmplam,0.1d0*alam)
    end do

    RETURN
    END SUBROUTINE lnsrch
