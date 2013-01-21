!---------------------------------------------------------!
!   Solve the equation using globally convergent          !
!   Newton-Raphson method                                 !
!---------------------------------------------------------!

    SUBROUTINE  SolveEquation(&
                &Nth, Nphi, Lmax, Lgrid, GLQWeights, GLQZeros,&
                &R, theta, phi,&
                &Maxit, tolA, tolF,&
                &LWORK,&
                &a)

        USE         omp_lib
        IMPLICIT    none
    
        !Declare passed variables
        INTEGER*4               Nth, Nphi, Lmax, Lgrid, Maxit
        INTEGER*4               LWORK
        REAL*8                  GLQWeights(Lgrid+1), GLQZeros(Lgrid+1)
        REAL*8                  R, theta(Nth), phi(Nphi)
        REAL*8                  tolA, tolF
        REAL*8                  a(2,Lmax+1,Lmax+1)

        !Declare local variables
        INTEGER*4               i, j, k, it
        REAL*8                  errF, errA

        REAL*8                  aV((Lmax+1)**2)
        REAL*8                  deltaA((Lmax+1)**2)
    
        REAL*8                  Res((Lmax+1)**2)

        REAL*8                  Jacobian((Lmax+1)**2,(Lmax+1)**2)
        REAL*8                  invJacobian((Lmax+1)**2,(Lmax+1)**2)

        !Needed for inverting the Jacobian
        REAL*8                  WORK(MAX(1,LWORK))
        INTEGER*4               IPIV((Lmax+1)**2)
        INTEGER*4               INFO

        REAL*8                  S((Lmax+1)**2)
        REAL*8                  U((Lmax+1)**2,(Lmax+1)**2)
        REAL*8                  VT((Lmax+1)**2,(Lmax+1)**2)

        !Main subroutine

        CALL SHCilmToVector(a,aV,Lmax)

        DO it=1,Maxit
        
            !Calculate the residual
            CALL RadiationFunction(Nth,Nphi,Lmax,Lgrid,&
                &GLQWeights,GLQZeros,R,theta,phi,&
                &aV,Res)

            !Compute the Jacobian by numerical differentiation
            CALL FDJacobian(Nth,Nphi,Lmax,Lgrid,&
                &GLQWeights,GLQZeros,R,theta,phi,&
                &aV,Jacobian)

            !Check the maximum value of the Jacobian vs. the minimum value
            PRINT *, 'Maximum value = ', MAXVAL(Jacobian)
            PRINT *, 'Minimum value = ', MINVAL(Jacobian)
            PRINT *, 'Ratio of max/min = ', MAXVAL(Jacobian)/MINVAL(Jacobian)

!!$            !Invert the Jacobian
!!$            invJacobian = Jacobian
!!$            CALL DGETRF( (Lmax+1)**2, (Lmax+1)**2, invJacobian,&
!!$                        &(Lmax+1)**2, IPIV, INFO)
!!$            IF( INFO .NE. 0 ) THEN
!!$                PRINT *, 'Error in LU factorization:'
!!$                PRINT *, 'Info =', INFO
!!$                STOP
!!$            END IF
!!$
!!$            CALL DGETRI( (Lmax+1)**2, invJacobian, (Lmax+1)**2,&
!!$                        &IPIV, WORK, LWORK, INFO)
!!$
!!$            IF( INFO .NE. 0 ) THEN
!!$                PRINT *, 'Error in inversion:'
!!$                PRINT *, 'Info =', INFO
!!$                STOP
!!$            END IF
!!$
!!$            IF( LWORK .EQ. -1 ) THEN
!!$                PRINT *, 'Work space query, optimal LWORK =', WORK(1)
!!$                STOP
!!$            END IF


            !Invert the Jacobian using Moore-Penrose inverse using SVD
            CALL DGESVD('S','S',(Lmax+1)**2,(Lmax+1)**2,Jacobian,&
                    &(Lmax+1)**2,S,U,(Lmax+1)**2,VT,(Lmax+1)**2,WORK,LWORK,INFO)

            IF( LWORK .EQ. -1 ) THEN
                PRINT *, 'Work space query, optimal LWORK =', WORK(1)
                STOP
            END IF

            !Compute invJ = VT**t * SIGMA * U**t
            !$OMP PARALLEL DO
            DO i=1,(Lmax+1)**2
                CALL DSCAL((Lmax+1)**2,(1.0D0/S(i)),U(1,i),1)
            END DO
            !$OMP END PARALLEL DO

            CALL DGEMM('C','C',(Lmax+1)**2,(Lmax+1)**2,(Lmax+1)**2,1.0D0,&
                &VT,(Lmax+1)**2,U,(Lmax+1)**2,0.0D0,invJacobian,(Lmax+1)**2)

            !Check if a is within tolerance
            deltaA = MATMUL(invJacobian,Res)
            errA = SQRT( DOT_PRODUCT(deltaA, deltaA) ) / ((Lmax+1)**2)
!!$            IF( errA .LE. tolA ) THEN
!!$                PRINT *, 'error in solution is within tolerance' 
!!$                CALL SHVectorToCilm(aV,a,Lmax)
!!$                RETURN
!!$            END IF

            !Update the solution
            aV = aV - deltaA

            !Check if residual is small enough
            errF = SQRT( DOT_PRODUCT(Res, Res) ) / ((Lmax+1)**2)
            IF( errF .LE. tolF ) THEN
                PRINT *, 'error in residual is within tolerance'
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
