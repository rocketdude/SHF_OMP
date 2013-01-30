!----------------------------------------------------------!
!     GetSchwarzschildMetric  subroutine                   !
!----------------------------------------------------------!

      SUBROUTINE GetSchwarzschildMetric(&
&Mass,&
&Nr, Nth, Nphi,&
&rho, theta, phi,&
&rmaxX, rmaxY, rmaxZ,&
&rminX, rminY, rminZ,&
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
        REAL*8               :: rho(Nr)
        REAL*8               :: theta(Nth)
        REAL*8               :: phi(Nphi)

        REAL*8               :: rmaxX,rmaxY,rmaxZ
        REAL*8               :: rminX,rminY,rminZ

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

        REAL*8                  alphaXYZ
        REAL*8                  betaX
        REAL*8                  betaY
        REAL*8                  betaZ
        REAL*8                  gXX
        REAL*8                  gYY
        REAL*8                  gZZ
        REAL*8                  gXY
        REAL*8                  gXZ
        REAL*8                  gYZ

        REAL*8                  gsph(4,4), gcart(4,4), JMatrix(4,4)

        REAL*8                    rmax, rmin
        REAL*8                    drmaxdth,drmindth,drmaxdphi,drmindphi
        REAL*8                    r(Nr)
        REAL*8                    x,y,z
        INTEGER*4                 i, j, k
        INTEGER*4                 error

!----------------------------------------------------------!
!      Main                                                !
!----------------------------------------------------------!

        PRINT *, 'Evaluating Metric with Mass = ', Mass

        !$OMP PARALLEL DO &
        !$OMP &PRIVATE(i,k,x,y,z,r,rmax,rmin,alphaXYZ,betaX,betaY,betaZ,&
        !$OMP &        gXX,gYY,gZZ,gXY,gXZ,gYZ,gsph,gcart,JMatrix,error,&
        !$OMP &        drmaxdth,drmindth,drmaxdphi,drmindphi)
        DO j = 1, Nth
            DO k = 1, Nphi

                 CALL EvaluateRadialExtent(rmaxX,rmaxY,rmaxZ,&
                                                theta(j),phi(k),rmax)
                 CALL EvaluateRadialExtent(rminX,rminY,rminZ,&
                                                theta(j),phi(k),rmin)
                 CALL GetRadialCoordinates(Nr,rmax,rmin,rho,r)

                IF( rmax .EQ. 0.0D0 ) THEN
                    drmaxdth = 0.0D0
                    drmaxdphi = 0.0D0
                ELSE
                    drmaxdth = ( ( rmaxX*COS(phi(k)) )**2 +&
                                &( rmaxY*SIN(phi(k)) )**2 -&
                                &  rmaxZ**2 ) * &
                                & SIN(theta(j))*COS(theta(j)) / rmax
                    drmaxdphi = ( rmaxY**2 - rmaxX**2 )*SIN(theta(j))**2 *&
                                & SIN(phi(k))*COS(phi(k)) / rmax
                END IF

                IF( rmin .EQ. 0.0D0 ) THEN
                    drmindth = 0.0D0
                    drmindphi = 0.0D0
                ELSE
                    drmindth = ( ( rminX*COS(phi(k)) )**2 +&
                                &( rminY*SIN(phi(k)) )**2 -&
                                &  rminZ**2 ) * &
                                & SIN(theta(j))*COS(theta(j)) / rmin
                    drmindphi = ( rminY**2 - rminX**2 )*SIN(theta(j))**2 *&
                                & SIN(phi(k))*COS(phi(k)) / rmin
                END IF

                DO i = 1, Nr

                    x = r(i)*SIN(theta(j))*COS(phi(k))
                    y = r(i)*SIN(theta(j))*SIN(phi(k))
                    z = r(i)*COS(theta(j))
                 
                    !Build the metric in cartesian coordinates
                    !g_ab = eta_ab + 2M/r l_a l_b
                    alphaXYZ = 1.0D0 / SQRT(1.0D0+2.0D0*Mass/r(i))
                    betaX    = (2.0D0*Mass*r(i)*x)/(r(i)**3+2.0D0*Mass*(r(i))**2)
                    betaY    = (2.0D0*Mass*r(i)*y)/(r(i)**3+2.0D0*Mass*(r(i))**2)
                    betaZ    = (2.0D0*Mass*r(i)*z)/(r(i)**3+2.0D0*Mass*(r(i))**2)
                    gXX      = 1.0D0 + 2.0D0*Mass*x*x/(r(i)**3)
                    gYY      = 1.0D0 + 2.0D0*Mass*y*y/(r(i)**3)
                    gZZ      = 1.0D0 + 2.0D0*Mass*z*z/(r(i)**3)
                    gXY      = 2.0D0*Mass*x*y/(r(i)**3)
                    gXZ      = 2.0D0*Mass*x*z/(r(i)**3)
                    gYZ      = 2.0D0*Mass*y*z/(r(i)**3)

                    !Build the Jacobian matrix JMatrix
                    CALL EvaluateJacobian(rmax,rmin,&
                                        &drmaxdth,drmindth,drmaxdphi,drmindphi,&
                                        &rho(i),theta(j),phi(k),&
                                        &JMatrix)

                    ! Build the matrix g_cart
                    CALL EvaluateMatrixofMetric(&
                         &alphaXYZ,&
                         &betaX, betaY, betaZ,&
                         &gXX, gYY, gZZ,&
                         &gXY, gXZ, gYZ,&
                         &gcart)

                    gsph = MATMUL(TRANSPOSE(JMatrix), MATMUL(gcart,JMatrix))

                    CALL Invert3Metric(gsph(2,2), gsph(3,3), gsph(4,4),&
                                    &gsph(2,3), gsph(2,4), gsph(3,4),&
                                    &gRR(i,j,k), gThTh(i,j,k), gPhiPhi(i,j,k),&
                                    &gRTh(i,j,k), gRPhi(i,j,k), gThPhi(i,j,k),&
                                    &error)

                    IF( error .NE. 0 ) THEN
                        PRINT *, 'Determinant is zero, error in metric invert'
                        PRINT *, 'i=',i,',j=',j,',k=',k
                        PRINT *, 'r=',r(i),',th=',theta(j),',phi=',phi(k)
                        PRINT *, 'gsph = ', gsph
                        STOP
                    END IF

                    ! Get beta^r, beta^th, beta^phi
                    ! gsph(1,2) = g_01 = beta_r, gsph(1,3) = g_02 = beta_th, 
                    ! gsph(1,4) = g_03 = beta_phi
                    betaR(i,j,k) = gRR(i,j,k)*gsph(1,2) &
                                &+ gRTh(i,j,k)*gsph(1,3) &
                                &+ gRPhi(i,j,k)*gsph(1,4)
                    betaTh(i,j,k)= gRTh(i,j,k)*gsph(1,2) &
                                &+ gThTh(i,j,k)*gsph(1,3) &
                                &+ gThPhi(i,j,k)*gsph(1,4)
                    betaPhi(i,j,k)= gRPhi(i,j,k)*gsph(1,2) &
                                &+ gThPhi(i,j,k)*gsph(1,3) &
                                &+ gPhiPhi(i,j,k)*gsph(1,4)

                    ! Get alpha
                    ! alpha = sqrt( beta^i beta_i - g_00 )
                    ! beta^r = betaR, beta^th = betaTh, beta^phi = betaPhi, 
                    ! gsph(1,1) = g_00
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
      END SUBROUTINE GetSchwarzschildMetric
