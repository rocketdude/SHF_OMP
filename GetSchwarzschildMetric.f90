!----------------------------------------------------------!
!     GetSchwarzschildMetric  subroutine                   !
!----------------------------------------------------------!

      SUBROUTINE GetSchwarzschildMetric(&
&Mass,&
&Nr, Nth, Nphi,&
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

        INTEGER*4            :: Nr,Nth,Nphi

        REAL*8               :: Mass
        REAL*8               :: r(Nr)
        REAL*8               :: theta(Nth)
        REAL*8               :: phi(Nphi)

        REAL*8, INTENT(out)  :: alpha(Nr,Nth,Nphi)
        REAL*8, INTENT(out)  :: betaR(Nr,Nth,Nphi)
        REAL*8, INTENT(out)  :: betaTh(Nr,Nth,Nphi)
        REAL*8, INTENT(out)  :: betaPhi(Nr,Nth,Nphi)

        REAL*8, INTENT(out)  :: gRR(Nr,Nth,Nphi)
        REAL*8, INTENT(out)  :: gThTh(Nr,Nth,Nphi)
        REAL*8, INTENT(out)  :: gPhiPhi(Nr,Nth,Nphi)

        REAL*8, INTENT(out)  :: gRTh(Nr,Nth,Nphi)
        REAL*8, INTENT(out)  :: gRPhi(Nr,Nth,Nphi)
        REAL*8, INTENT(out)  :: gThPhi(Nr,Nth,Nphi)

!-----------------------------------------------------------!
!       Declare local variables                             !
!-----------------------------------------------------------!

        INTEGER*4                 i, j, k

!----------------------------------------------------------!
!      Main                                                !
!----------------------------------------------------------!

        PRINT *, 'Evaluating Metric'

        !$OMP PARALLEL DO SHARED(Nr, Nth, Nphi, Mass, r, theta, phi)&
        !$OMP &PRIVATE(j, k)
        DO i = 1, Nr
           DO j = 1, Nth
              DO k = 1, Nphi

                 alpha(i,j,k) = SQRT( r(i) / ( r(i) + 2.0D0*Mass ) )
                 
                 betaR(i,j,k) = 2.0D0*Mass / (2.0D0*Mass + r(i) )
                 betaTh(i,j,k) = 0.0D0
                 betaPhi(i,j,k) = 0.0D0

                 gRR(i,j,k) = 1.0D0 / ( 1.0D0 + 2.0D0*Mass / r(i) )
                 gThTh(i,j,k) = 1.0D0/(r(i)**2)
                 gPhiPhi(i,j,k) = 1.0D0/( r(i) * SIN( theta(j) ) )**2

                 gRTh(i,j,k) = 0.0D0
                 gRPhi(i,j,k) = 0.0D0
                 gThPhi(i,j,k) = 0.0D0

              END DO
           END DO
        END DO
        !$OMP END PARALLEL DO

        PRINT *, 'DONE!'

        RETURN
      END SUBROUTINE GetSchwarzschildMetric
