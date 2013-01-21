!---------------------------------------------------------!
!   Solve the equation using Newton-Raphson               !
!---------------------------------------------------------!

    SUBROUTINE  SolveEquation(&
                &Nth, Nphi, Lmax, Lgrid, GLQWeights, GLQZeros,&
                &R, theta, phi,&
                &Maxit, tolA, tolF,&
                &p, eps,&
                &LWORK,&
                &a)

        USE         omp_lib
        IMPLICIT    none
    
        !Declare passed variables
        INTEGER*4               Nth, Nphi, Lmax, Lgrid, Maxit
        INTEGER*4               p
        INTEGER*4               LWORK
        REAL*8                  GLQWeights(Lgrid+1), GLQZeros(Lgrid+1)
        REAL*8                  R, theta(Nth), phi(Nphi)
        REAL*8                  tolA, tolF
        REAL*8                  eps
        REAL*8                  a(2,Lmax+1,Lmax+1)

        !Declare local variables
        INTEGER*4               YilmIndex
        INTEGER*4               i, j, k, it, l, ml, ccol
        REAL*8                  MaxDA
        REAL*8                  errF, errA
 
        REAL*8                  Filter(Lmax+1)
        
        REAL*8                  aVector((Lmax+1)**2)
        REAL*8                  deltaA((Lmax+1)**2)
    
        REAL*8                  Res(2,Lmax+1,Lmax+1)
        REAL*8                  ResVector((Lmax+1)**2)

        !Needed for calculating the Jacobian
        REAL*8                  aEval(2,Lmax+1,Lmax+1,4)

        REAL*8                  da(2,Lmax+1,Lmax+1)
        REAL*8                  i12da(2,Lmax+1,Lmax+1) 

        REAL*8                  Fn(2,Lmax+1,Lmax+1,4)
        REAL*8                  dFnda(2,Lmax+1,Lmax+1)
        REAL*8                  dFndaVector((Lmax+1)**2)
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

        CALL GetFilter(Lmax,p,eps,Filter)

        MaxDA = 1.0D-3

        DO it=1,Maxit
        
            !Compute the Jacobian by numerical differentiation
            !$OMP PARALLEL DO PRIVATE(i)
            DO j=1,(Lmax+1)
                DO i=1,(Lmax+1)
                    da(1,i,j) = MIN( MaxDA*a(1,i,j), MaxDA )
                    da(2,i,j) = MIN( MaxDA*a(2,i,j), MaxDA )

                    IF( da(1,i,j) .NE. 0.0D0 ) THEN
                        i12da(1,i,j) = 1.0D0/(12.0D0*da(1,i,j))
                    ELSE
                        i12da(1,i,j) = 0.0D0
                    END IF

                    IF( da(2,i,j) .NE. 0.0D0 ) THEN
                        i12da(2,i,j) = 1.0D0/(12.0D0*da(2,i,j))
                    ELSE
                        i12da(2,i,j) = 0.0D0
                    END IF
                END DO
            END DO
            !$OMP END PARALLEL DO

            aEval(:,:,:,1) = a
            aEval(:,:,:,2) = a
            aEval(:,:,:,3) = a
            aEval(:,:,:,4) = a

            !Not threadsafe
            DO l=0,Lmax

                i = 1
                DO ml=0,l

                    !Setup differentiation grid
                    aEval(i,l+1,ml+1,1) = a(i,l+1,ml+1) - 2.0D0*da(i,l+1,ml+1)
                    aEval(i,l+1,ml+1,2) = a(i,l+1,ml+1) - da(i,l+1,ml+1)
                    aEval(i,l+1,ml+1,3) = a(i,l+1,ml+1) + da(i,l+1,ml+1)
                    aEval(i,l+1,ml+1,4) = a(i,l+1,ml+1) + 2.0D0*da(i,l+1,ml+1)

                    !Evaluate the function at the different "grid" points
                    CALL RadiationFunction(Nth,Nphi,Lmax,Lgrid,&
                        &GLQWeights,GLQZeros,R,theta,phi,&
                        &aEval(:,:,:,1),Fn(:,:,:,1))

                    CALL RadiationFunction(Nth,Nphi,Lmax,Lgrid,&
                        &GLQWeights,GLQZeros,R,theta,phi,&
                        &aEval(:,:,:,2),Fn(:,:,:,2))

                    CALL RadiationFunction(Nth,Nphi,Lmax,Lgrid,&
                        &GLQWeights,GLQZeros,R,theta,phi,&
                        &aEval(:,:,:,3),Fn(:,:,:,3))

                    CALL RadiationFunction(Nth,Nphi,Lmax,Lgrid,&
                        &GLQWeights,GLQZeros,R,theta,phi,&
                        &aEval(:,:,:,4),Fn(:,:,:,4))

                    !Calculate numerical derivatives
                    dFnda = i12da*&
                        &( Fn(:,:,:,1) - 8.0D0*Fn(:,:,:,2) +&
                           8.0D0*Fn(:,:,:,3) - Fn(:,:,:,4) )

                    !Transfer to Jacobian
                    CALL SHCilmToVector(dFnda,dFndaVector,Lmax)
                    ccol = YilmIndex(i,l,ml)
                    Jacobian(:,ccol) = dFndaVector

                    !Restore aEval
                    aEval(i,l+1,ml+1,:) = a(i,l+1,ml+1)

                END DO

                i = 2
                DO ml=1,l

                    !Setup differentiation grid
                    aEval(i,l+1,ml+1,1) = a(i,l+1,ml+1) - 2.0D0*da(i,l+1,ml+1)
                    aEval(i,l+1,ml+1,2) = a(i,l+1,ml+1) - da(i,l+1,ml+1)
                    aEval(i,l+1,ml+1,3) = a(i,l+1,ml+1) + da(i,l+1,ml+1)
                    aEval(i,l+1,ml+1,4) = a(i,l+1,ml+1) + 2.0D0*da(i,l+1,ml+1)

                    !Evaluate the function at the different "grid" points
                    CALL RadiationFunction(Nth,Nphi,Lmax,Lgrid,&
                        &GLQWeights,GLQZeros,R,theta,phi,&
                        &aEval(:,:,:,1),Fn(:,:,:,1))

                    CALL RadiationFunction(Nth,Nphi,Lmax,Lgrid,&
                        &GLQWeights,GLQZeros,R,theta,phi,&
                        &aEval(:,:,:,2),Fn(:,:,:,2))

                    CALL RadiationFunction(Nth,Nphi,Lmax,Lgrid,&
                        &GLQWeights,GLQZeros,R,theta,phi,&
                        &aEval(:,:,:,3),Fn(:,:,:,3))

                    CALL RadiationFunction(Nth,Nphi,Lmax,Lgrid,&
                        &GLQWeights,GLQZeros,R,theta,phi,&
                        &aEval(:,:,:,4),Fn(:,:,:,4))

                    !Calculate numerical derivatives
                    dFnda = i12da*&
                        &( Fn(:,:,:,1) - 8.0D0*Fn(:,:,:,2) +&
                           8.0D0*Fn(:,:,:,3) - Fn(:,:,:,4) )

                    !Transfer to Jacobian
                    CALL SHCilmToVector(dFnda,dFndaVector,Lmax)
                    ccol = YilmIndex(i,l,ml)
                    Jacobian(:,ccol) = dFndaVector

                    !Restore aEval
                    aEval(i,l+1,ml+1,:) = a(i,l+1,ml+1)

                END DO
            END DO

            !Calculate the residual
            CALL RadiationFunction(Nth,Nphi,Lmax,Lgrid,&
                &GLQWeights,GLQZeros,R,theta,phi,&
                &a,Res)
            CALL SHCilmToVector(Res,ResVector,Lmax)

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
            deltaA = MATMUL(invJacobian,ResVector)
            errA = SQRT( DOT_PRODUCT(deltaA, deltaA) ) / ((Lmax+1)**2)
            !IF( errA .LE. tolA ) THEN
            !    PRINT *, 'error in solution is within tolerance' 
            !    RETURN
            !END IF

            !Update the solution
            CALL SHCilmToVector(a,aVector,Lmax)
            aVector = aVector - deltaA
            CALL SHVectorToCilm(aVector,a,Lmax)

            !Check if residual is small enough
            errF = SQRT( DOT_PRODUCT(ResVector, ResVector) ) / ((Lmax+1)**2)
            IF( errF .LE. tolF ) THEN
                PRINT *, 'error in residual is within tolerance'
                RETURN
            END IF

            !$OMP PARALLEL DO
            DO l=0, Lmax
                a(:,l+1,:) = a(:,l+1,:) * Filter(l+1)
            END DO
            !$OMP END PARALLEL DO

            !Output to user
            PRINT *, 'Iteration # =', it
            PRINT *, 'error in residual =', errF
            PRINT *, 'error in solution =', errA
            PRINT *, '========================='
            
        END DO

        PRINT *, 'ERROR: Maximum number of iterations reached'
        RETURN
    END SUBROUTINE
