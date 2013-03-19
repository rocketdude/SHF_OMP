!----------------------------------------------------------!
!     GetSchwarzschildMetric  subroutine                   !
!----------------------------------------------------------!

      SUBROUTINE GetKerrSchildMetric(&
&Mass, Spin,&
&Nr, Nth, Nphi,&
&rmax, rmin,&
&rho, theta, phi,&
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

        REAL*8               :: Mass, Spin
        REAL*8               :: rmax, rmin
        REAL*8               :: rho(Nr)
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
        REAL*8                    r(Nr)
        REAL*8                    H,rBL,x,y,z
        REAL*8                    JMatrix(4,4),gcart(4,4),gsph(4,4)
        REAL*8                    l1,l2,l3
        REAL*8                    alp,bx,by,bz,gxx,gyy,gzz,gxy,gxz,gyz
        INTEGER*4                 error

!----------------------------------------------------------!
!      Main                                                !
!----------------------------------------------------------!

        PRINT *, 'Evaluating Metric'
        CALL GetRadialCoordinates(Nr,rmax,rmin,rho,r)

        !$OMP PARALLEL DO SHARED(Nr, Nth, Nphi, Mass, Spin, r, theta, phi)&
        !$OMP &PRIVATE(j,k,l1,l2,l3,gxx,gyy,gzz,gxy,gxz,gyz,alp,bx,by,bz,&
        !$OMP & JMatrix,gcart,gsph,H,rBL,x,y,z,error)
        DO i = 1, Nr
           DO j = 1, Nth
              DO k = 1, Nphi
                
                x = r(i)*SIN(theta(j))*COS(phi(k))
                y = r(i)*SIN(theta(j))*SIN(phi(k))
                z = r(i)*COS(theta(j))

                rBL = SQRT( r(i)*r(i)-Spin*Spin +&
                        & SQRT((2.0D0*Spin*z)**2 + (Spin*Spin-r(i)*r(i))**2))&
                        / SQRT(2.0D0)

                H = Mass*(rBL**3)/(rBL**4 + Spin*Spin*z*z)
                l1 = (rBL*x+Spin*y)/(rBL*rBL+Spin*Spin)
                l2 = (rBL*y-Spin*x)/(rBL*rBL+Spin*Spin)
                l3 = z/rBL

                gxx = 1.0D0 + 2.0D0*H*l1*l1
                gyy = 1.0D0 + 2.0D0*H*l2*l2
                gzz = 1.0D0 + 2.0D0*H*l3*l3
                gxy = 2.0D0*H*l1*l2
                gxz = 2.0D0*H*l1*l3
                gyz = 2.0D0*H*l2*l3

                alp = 1.0D0 / SQRT( 1.0D0 + 2.0D0*H )
                bx  = 2.0D0*H*l1 / (1.0D0 + 2.0D0*H)
                by  = 2.0D0*H*l2 / (1.0D0 + 2.0D0*H)
                bz  = 2.0D0*H*l3 / (1.0D0 + 2.0D0*H)

                CALL EvaluateJacobian(rmax,rmin,rho(i),theta(j),phi(k),JMatrix)
                CALL EvaluateMatrixofMetric(alp,bx,by,bz,&
                                           &gxx,gyy,gzz,gxy,gxz,gyz,gcart)
                gsph = MATMUL(TRANSPOSE(JMatrix),MATMUL(gcart,JMatrix))
                
                CALL Invert3Metric(gsph(2,2), gsph(3,3), gsph(4,4),&
                                    &gsph(2,3), gsph(2,4), gsph(3,4),&
                                    &gRR(i,j,k), gThTh(i,j,k), gPhiPhi(i,j,k),&
                                    &gRTh(i,j,k), gRPhi(i,j,k), gThPhi(i,j,k),&
                                    &error)

                IF(error.NE.0) STOP "***ERROR in metric inversion***"

                    betaR(i,j,k) = gRR(i,j,k)*gsph(1,2) &
                                &+ gRTh(i,j,k)*gsph(1,3) &
                                &+ gRPhi(i,j,k)*gsph(1,4)
                    betaTh(i,j,k)= gRTh(i,j,k)*gsph(1,2) &
                                &+ gThTh(i,j,k)*gsph(1,3) &
                                &+ gThPhi(i,j,k)*gsph(1,4)
                    betaPhi(i,j,k)= gRPhi(i,j,k)*gsph(1,2) &
                                &+ gThPhi(i,j,k)*gsph(1,3) &
                                &+ gPhiPhi(i,j,k)*gsph(1,4)
                    alpha(i,j,k) = SQRT( gsph(1,2) * betaR(i,j,k) + &
                                    &gsph(1,3) * betaTh(i,j,k) + &
                                    &gsph(1,4) * betaPhi(i,j,k) &
                                & - gsph(1,1) )
                

              END DO
           END DO
        END DO
        !$OMP END PARALLEL DO

        PRINT *, 'DONE!'

        RETURN
      END SUBROUTINE GetKerrSchildMetric
