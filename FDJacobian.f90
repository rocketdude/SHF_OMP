!---------------------------------------------------------!
!   Compute the Jacobian using Finite-Difference          !
!---------------------------------------------------------!

    SUBROUTINE FDJacobian(Nth,Nphi,Lmax,Lgrid,&
                        & GLQWeights,GLQZeros,R,theta,phi,&
                        & a,Jacobian)

        !Computes the Jacobian using central difference of O(h^4)

        USE         omp_lib
        IMPLICIT    none

        !Declare passed variables
        INTEGER*4               Nth,Nphi,Lmax,Lgrid
        INTEGER*4               GLQWeights(Lgrid+1), GLQZeros(Lgrid+1)
        REAL*8                  R, theta(Nth), phi(Nphi)
        REAL*8                  a((Lmax+1)**2)
        REAL*8,INTENT(OUT)  ::  Jacobian((Lmax+1)**2,(Lmax+1)**2)

        !Declare Local variables
        REAL*8                  MaxDA
        REAL*8                  aEval((Lmax+1)**2,4)
        REAL*8                  da((Lmax+1)**2)
        REAL*8                  i12da((Lmax+1)**2)
        REAL*8                  Fn((Lmax+1)**2,4)
        REAL*8                  dFnda((Lmax+1)**2)
        INTEGER*4               i,j,k,N

        !Main subroutine

        N = (Lmax+1)**2
        MaxDA = 1.0D-3

        !$OMP PARALLEL DO
        DO i=1,N
            da(i) = MIN(MaxDA*a(i), MaxDA)

            IF( da(i) .NE. 0.0D0 ) THEN
                i12da(i) = 1.0D0/(12.0D0*da(i))
            ELSE
                i12da(i) = 0.0D0
            END IF
        END DO
        !$OMP END PARALLEL DO

        aEval(:,1) = a
        aEval(:,2) = a
        aEval(:,3) = a
        aEval(:,4) = a

        !This is not threadsafe
        DO i=1,N
            !Setup differentiation grid
            aEval(i,1) = a(i) - 2.0D0*da(i)
            aEval(i,2) = a(i) - da(i)
            aEval(i,3) = a(i) + da(i)
            aEval(i,4) = a(i) + 2.0D0*da(i)

            !Evaluate the function at different "grid" points
            CALL RadiationFunction(Nth,Nphi,Lmax,Lgrid,GLQWeights,GLQZeros,&
                                    &R,theta,phi,&
                                    &aEval(:,1),Fn(:,1))
            CALL RadiationFunction(Nth,Nphi,Lmax,Lgrid,GLQWeights,GLQZeros,&
                                    &R,theta,phi,&
                                    &aEval(:,2),Fn(:,2))
            CALL RadiationFunction(Nth,Nphi,Lmax,Lgrid,GLQWeights,GLQZeros,&
                                    &R,theta,phi,&
                                    &aEval(:,3),Fn(:,3))
            CALL RadiationFunction(Nth,Nphi,Lmax,Lgrid,GLQWeights,GLQZeros,&
                                    &R,theta,phi,&
                                    &aEval(:,4),Fn(:,4))

            !Calculate numerical derivatives
            dFnda = i12da*(Fn(:,1) - 8.0D0*Fn(:,2) + 8.0D0*Fn(:,3) - Fn(:,4))

            !Transfer to Jacobian
            Jacobian(:,i) = dFnda

            !Restore aEval
            aEval(i,:) = a(i)
        END DO

    RETURN
    END SUBROUTINE FDJacobian

