!----------------------------------------------------------!
!     GetSchwarzschildMetric  subroutine                   !
!----------------------------------------------------------!

      SUBROUTINE GetSchwarzschildMetric(&
&Mass,&
&M, Mr, NP,&
&r, theta, phi,&
&alpha,&
&betaR, betaTh, betaPhi,&
&gRR, gThTh, gPhiPhi,&
&gRTh, gRPhi, gThPhi)

        !This subroutine calculates the Schwarzschild metric
        !with all indices up (for both g's and beta's)
        !For Minkowski metric, put BHoleMass = 0

        USE       omp_lib
        IMPLICIT  none

!-----------------------------------------------------------!
!       Declare calling variables                           !
!-----------------------------------------------------------!

        INTEGER*4            :: M, Mr, NP

        REAL*8               :: Mass
        REAL*8               :: r(Mr+1)
        REAL*8               :: theta(2*M)
        REAL*8               :: phi(2*M)

        REAL*8, INTENT(out)  :: alpha(4*NP)
        REAL*8, INTENT(out)  :: betaR(4*NP)
        REAL*8, INTENT(out)  :: betaTh(4*NP)
        REAL*8, INTENT(out)  :: betaPhi(4*NP)

        REAL*8, INTENT(out)  :: gRR(4*NP)
        REAL*8, INTENT(out)  :: gThTh(4*NP)
        REAL*8, INTENT(out)  :: gPhiPhi(4*NP)

        REAL*8, INTENT(out)  :: gRTh(4*NP)
        REAL*8, INTENT(out)  :: gRPhi(4*NP)
        REAL*8, INTENT(out)  :: gThPhi(4*NP)

!-----------------------------------------------------------!
!       Declare local variables                             !
!-----------------------------------------------------------!

        INTEGER*4                 i, j, k
        INTEGER*4                 crow    !Row counter

!----------------------------------------------------------!
!      Main                                                !
!----------------------------------------------------------!

        PRINT *, 'Evaluating Metric'

        !$OMP PARALLEL DO SHARED(M, Mass, r, theta, phi)&
        !$OMP &PRIVATE(crow, j, k)
        DO i = 1, Mr+1
           DO j = 1, 2*M
              DO k = 1, 2*M

                 crow = (i-1)*(2*M)*(2*M) + (j-1)*(2*M) + k

                 alpha(crow) = SQRT( r(i) / ( r(i) + 2.0D0*Mass ) )
                 
                 betaR(crow) = 2.0D0*Mass / (2.0D0*Mass + r(i) )
                 betaTh(crow) = 0.0D0
                 betaPhi(crow) = 0.0D0

                 gRR(crow) = 1.0D0 / ( 1.0D0 + 2.0D0*Mass / r(i) )
                 gThTh(crow) = 1.0D0/(r(i)**2)
                 gPhiPhi(crow) = 1.0D0/( r(i) * SIN( theta(j) ) )**2

                 gRTh(crow) = 0.0D0
                 gRPhi(crow) = 0.0D0
                 gThPhi(crow) = 0.0D0

              END DO
           END DO
        END DO
        !$OMP END PARALLEL DO

        PRINT *, 'DONE!'

        RETURN
      END SUBROUTINE GetSchwarzschildMetric
