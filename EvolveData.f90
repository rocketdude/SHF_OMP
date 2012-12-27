!--------------------------------------------------------!
!    EvolveData subroutine                               !
!--------------------------------------------------------!

      SUBROUTINE EvolveData(&
&Nr, Nth, Nphi, Mr, Lmax, Lgrid,&
&GLQWeights, GLQZeros,&
&rmin, rmax,&
&r, rho, theta, phi,&
&epsStar, alpStar,&
&SolPhi,&
&t, dt,&
&a)

        !EvolveData subroutine uses the SSPRK (5,4)
        !An outline of the integrator can be found on 
        !"Spectral Methods for Time-Dependent Problems" 
        !by Jan Hesthaven, pg 201
        !with CFL coefficient of 1.508
        
        USE             omp_lib
        IMPLICIT        none

!--------------------------------------------------------!
!     Declare calling variables                          !
!--------------------------------------------------------!

        INTEGER*4               :: Nr,Nth,Nphi,Mr,Lmax,Lgrid

        REAL*8                  :: rmin
        REAL*8                  :: rmax
        REAL*8                  :: GLQWeights(Lgrid+1), GLQZeros(Lgrid+1)
        REAL*8                  :: r(Nr), rho(Nr), theta(Nth), phi(Nphi)

        REAL*8                  :: epsStar, alpStar
        REAL*8                  :: SolPhi

        REAL*8                  :: t
        REAL*8                  :: dt

        COMPLEX*16              :: a(Mr+1,2,Lmax+1,Lmax+1)

!--------------------------------------------------------!
!     Declare Locals                                     !
!--------------------------------------------------------!

        CHARACTER*32    CTemp
        INTEGER*4         i, j, k
        INTEGER*4         n
        
        REAL*8          PI

        REAL*8          rsqrd
        REAL*8          bdryterm
        REAL*8          sigma           !Stefan-Boltzmann constant

        COMPLEX*16      S(Nr,Nth,Nphi)

        !First spatial derivatives
        COMPLEX*16      dSdr(Nr,Nth,Nphi)
        COMPLEX*16      dSdth(Nr,Nth,Nphi)
        COMPLEX*16      dSdphi(Nr,Nth,Nphi)
        !Second spatical derivatives
        COMPLEX*16      d2Sdr2(Nr,Nth,Nphi)
        COMPLEX*16      d2Sdth2(Nr,Nth,Nphi)
        COMPLEX*16      d2Sdphi2(Nr,Nth,Nphi)

        COMPLEX*16      dSdt(Nr,Nth,Nphi)

        COMPLEX*16      aaRK(Mr+1,2,Lmax+1,Lmax+1,5)
        COMPLEX*16      dadt(Mr+1,2,Lmax+1,Lmax+1)
        COMPLEX*16      dadt3(Mr+1,2,Lmax+1,Lmax+1)

!--------------------------------------------------------!
!      Main Subroutine                                   !
!--------------------------------------------------------!

        PI = 3.141592653589793238462643383279502884197D0
        sigma = 5.670373D-8

!--------------------------------------------------------!
!      STEP 1                                            !
!--------------------------------------------------------!

        !Store the coefficients inside aaRK 
        !which will be used in the Runge-Kutta method
        !$OMP PARALLEL DO
        DO n=1,(Mr+1)
            aaRK(n,:,:,:,1) = a(n,:,:,:)
        END DO
        !$OMP END PARALLEL DO

        CALL SpectralToSpatialTransform(Nr,Nth,Nphi,Mr,Lmax,Lgrid,&
                &GLQWeights,GLQZeros,rho,theta,phi,a,S)

        !Calculate the derivatives
        CALL EvaluatedSdr(Nr,Nth,Nphi,Mr,Lmax,Lgrid,GLQWeights,GLQZeros,&
                &rho,theta,phi,rmax,rmin,a,dSdr)
        CALL EvaluatedSdphi(Nr,Nth,Nphi,Mr,Lmax,Lgrid,GLQWeights,GLQZeros,&
                &rho,theta,phi,a,dSdphi)
        CALL EvaluatedSdtheta(Nr,Nth,Nphi,Mr,Lmax,Lgrid,GLQWeights,GLQZeros,&
                &rho,theta,phi,a,dSdth)
        CALL Evaluated2Sdr2(Nr,Nth,Nphi,Mr,Lmax,Lgrid,GLQWeights,GLQZeros,&
                &rho,theta,phi,rmax,rmin,dSdr,d2Sdr2)
        CALL Evaluated2Sdphi2(Nr,Nth,Nphi,Mr,Lmax,Lgrid,GLQWeights,GLQZeros,&
                &rho,theta,phi,dSdphi,d2Sdphi2)
        CALL Evaluated2Sdtheta2(Nr,Nth,Nphi,Mr,Lmax,Lgrid,GLQWeights,GLQZeros,&
                &rho,theta,phi,dSdth,d2Sdth2)
 
        !$OMP PARALLEL DO PRIVATE(j, k, rsqrd, bdryterm)
        DO i = 1,Nr
           DO j = 1,Nth
              DO k = 1,Nphi
                 
                 IF( i .EQ. 1 ) THEN
                    !BOUNDARIES (APPLY BOUNDARY CONDITIONS)
                    rsqrd = r(i) * r(i)
                    bdryterm = -1.0D0*epsStar*sigma*( S(i,j,k)**4 )

                    IF( theta(j) .LT. PI/2.0D0 ) THEN
                        bdryterm = bdryterm + &
                            &alpStar*SolPhi*COS(theta(j))
                    END IF

                    dSdt(i,j,k) = &
                        &2.0D0*bdryterm/r(i) - &
                        &4.0D0*epsStar*sigma*bdryterm*( S(i,j,k)**3 ) + &
                        &dSdth(i,j,k)/(rsqrd*TAN(theta(j))) + &
                        &d2Sdth2(i,j,k)/rsqrd + &
                        &d2Sdphi2(i,j,k)/(rsqrd*SIN(theta(j))*SIN(theta(j)))

                 ELSE
                    !INNER POINTS
                    rsqrd = r(i) * r(i)

                    dSdt(i,j,k) = &
                        &2.0D0*dSdr(i,j,k)/r(i) + d2Sdr2(i,j,k) + &
                        &dSdth(i,j,k)/(rsqrd*TAN(theta(j))) + &
                        &d2Sdth2(i,j,k)/rsqrd + &
                        &d2Sdphi2(i,j,k)/(rsqrd*SIN(theta(j))*SIN(theta(j)))

                 END IF

              END DO
           END DO
        END DO
        !$OMP END PARALLEL DO
        
        CALL SpatialToSpectralTransform(Nr,Nth,Nphi,Mr,Lmax,Lgrid,&
                                       &GLQWeights,GLQZeros,&
                                       &rho,theta,phi,&
                                       &dSdt,dadt)

        !$OMP PARALLEL DO
        DO n=1,(Mr+1)
            a(n,:,:,:) = aaRK(n,:,:,:,1) + 0.391752226571890D0 * dt &
                    &* dadt(n,:,:,:)
        END DO
        !$OMP END PARALLEL DO


!--------------------------------------------------------!
!      STEP 2                                            !
!--------------------------------------------------------!

        !Store the coefficients inside aaRK 
        !which will be used in the Runge-Kutta method
        !$OMP PARALLEL DO
        DO n=1,(Mr+1)
            aaRK(n,:,:,:,2) = a(n,:,:,:)
        END DO
        !$OMP END PARALLEL DO
 
        CALL SpectralToSpatialTransform(Nr,Nth,Nphi,Mr,Lmax,Lgrid,&
                &GLQWeights,GLQZeros,rho,theta,phi,a,S)

        !Calculate the derivatives
        CALL EvaluatedSdr(Nr,Nth,Nphi,Mr,Lmax,Lgrid,GLQWeights,GLQZeros,&
                &rho,theta,phi,rmax,rmin,a,dSdr)
        CALL EvaluatedSdphi(Nr,Nth,Nphi,Mr,Lmax,Lgrid,GLQWeights,GLQZeros,&
                &rho,theta,phi,a,dSdphi)
        CALL EvaluatedSdtheta(Nr,Nth,Nphi,Mr,Lmax,Lgrid,GLQWeights,GLQZeros,&
                &rho,theta,phi,a,dSdth)
        CALL Evaluated2Sdr2(Nr,Nth,Nphi,Mr,Lmax,Lgrid,GLQWeights,GLQZeros,&
                &rho,theta,phi,rmax,rmin,dSdr,d2Sdr2)
        CALL Evaluated2Sdphi2(Nr,Nth,Nphi,Mr,Lmax,Lgrid,GLQWeights,GLQZeros,&
                &rho,theta,phi,dSdphi,d2Sdphi2)
        CALL Evaluated2Sdtheta2(Nr,Nth,Nphi,Mr,Lmax,Lgrid,GLQWeights,GLQZeros,&
                &rho,theta,phi,dSdth,d2Sdth2)
        
        !$OMP PARALLEL DO PRIVATE(j, k, rsqrd, bdryterm)
        DO i = 1,Nr
           DO j = 1,Nth
              DO k = 1,Nphi
                 
                 IF( i .EQ. 1 ) THEN
                    !BOUNDARIES (APPLY BOUNDARY CONDITIONS)
                    rsqrd = r(i) * r(i)
                    bdryterm = -1.0D0*epsStar*sigma*( S(i,j,k)**4 )

                    IF( theta(j) .LT. PI/2.0D0 ) THEN
                        bdryterm = bdryterm + &
                            &alpStar*SolPhi*COS(theta(j))
                    END IF

                    dSdt(i,j,k) = &
                        &2.0D0*bdryterm/r(i) - &
                        &4.0D0*epsStar*sigma*bdryterm*( S(i,j,k)**3 ) + &
                        &dSdth(i,j,k)/(rsqrd*TAN(theta(j))) + &
                        &d2Sdth2(i,j,k)/rsqrd + &
                        &d2Sdphi2(i,j,k)/(rsqrd*SIN(theta(j))*SIN(theta(j)))

                 ELSE
                    !INNER POINTS
                    rsqrd = r(i) * r(i)

                    dSdt(i,j,k) = &
                        &2.0D0*dSdr(i,j,k)/r(i) + d2Sdr2(i,j,k) + &
                        &dSdth(i,j,k)/(rsqrd*TAN(theta(j))) + &
                        &d2Sdth2(i,j,k)/rsqrd + &
                        &d2Sdphi2(i,j,k)/(rsqrd*SIN(theta(j))*SIN(theta(j)))

                 END IF

              END DO
           END DO
        END DO
        !$OMP END PARALLEL DO

        CALL SpatialToSpectralTransform(Nr,Nth,Nphi,Mr,Lmax,Lgrid,&
                                       &GLQWeights,GLQZeros,&
                                       &rho,theta,phi,&
                                       &dSdt,dadt)

        !$OMP PARALLEL DO
        DO n=1,(Mr+1)
            a(n,:,:,:) = 0.444370493651235D0 * aaRK(n,:,:,:,1) + &
              & 0.555629506348765D0 * aaRK(n,:,:,:,2) + &
              & 0.368410593050371D0 * dt * dadt(n,:,:,:)
        END DO
        !$OMP END PARALLEL DO


!--------------------------------------------------------!
!      STEP 3                                            !
!--------------------------------------------------------!

        !Store the coefficients inside aaRK 
        !which will be used in the Runge-Kutta method
        !$OMP PARALLEL DO
        DO n=1,(Mr+1)
            aaRK(n,:,:,:,3) = a(n,:,:,:)
        END DO
        !$OMP END PARALLEL DO

        CALL SpectralToSpatialTransform(Nr,Nth,Nphi,Mr,Lmax,Lgrid,&
                &GLQWeights,GLQZeros,rho,theta,phi,a,S)

        !Calculate the derivatives
        CALL EvaluatedSdr(Nr,Nth,Nphi,Mr,Lmax,Lgrid,GLQWeights,GLQZeros,&
                &rho,theta,phi,rmax,rmin,a,dSdr)
        CALL EvaluatedSdphi(Nr,Nth,Nphi,Mr,Lmax,Lgrid,GLQWeights,GLQZeros,&
                &rho,theta,phi,a,dSdphi)
        CALL EvaluatedSdtheta(Nr,Nth,Nphi,Mr,Lmax,Lgrid,GLQWeights,GLQZeros,&
                &rho,theta,phi,a,dSdth)
        CALL Evaluated2Sdr2(Nr,Nth,Nphi,Mr,Lmax,Lgrid,GLQWeights,GLQZeros,&
                &rho,theta,phi,rmax,rmin,dSdr,d2Sdr2)
        CALL Evaluated2Sdphi2(Nr,Nth,Nphi,Mr,Lmax,Lgrid,GLQWeights,GLQZeros,&
                &rho,theta,phi,dSdphi,d2Sdphi2)
        CALL Evaluated2Sdtheta2(Nr,Nth,Nphi,Mr,Lmax,Lgrid,GLQWeights,GLQZeros,&
                &rho,theta,phi,dSdth,d2Sdth2)
        
        !$OMP PARALLEL DO PRIVATE(j, k, rsqrd, bdryterm)
        DO i = 1,Nr
           DO j = 1,Nth
              DO k = 1,Nphi
                 
                 IF( i .EQ. 1 ) THEN
                    !BOUNDARIES (APPLY BOUNDARY CONDITIONS)
                    rsqrd = r(i) * r(i)
                    bdryterm = -1.0D0*epsStar*sigma*( S(i,j,k)**4 )

                    IF( theta(j) .LT. PI/2.0D0 ) THEN
                        bdryterm = bdryterm + &
                            &alpStar*SolPhi*COS(theta(j))
                    END IF

                    dSdt(i,j,k) = &
                        &2.0D0*bdryterm/r(i) - &
                        &4.0D0*epsStar*sigma*bdryterm*( S(i,j,k)**3 ) + &
                        &dSdth(i,j,k)/(rsqrd*TAN(theta(j))) + &
                        &d2Sdth2(i,j,k)/rsqrd + &
                        &d2Sdphi2(i,j,k)/(rsqrd*SIN(theta(j))*SIN(theta(j)))

                 ELSE
                    !INNER POINTS
                    rsqrd = r(i) * r(i)

                    dSdt(i,j,k) = &
                        &2.0D0*dSdr(i,j,k)/r(i) + d2Sdr2(i,j,k) + &
                        &dSdth(i,j,k)/(rsqrd*TAN(theta(j))) + &
                        &d2Sdth2(i,j,k)/rsqrd + &
                        &d2Sdphi2(i,j,k)/(rsqrd*SIN(theta(j))*SIN(theta(j)))

                 END IF

              END DO
           END DO
        END DO
        !$OMP END PARALLEL DO
 
        CALL SpatialToSpectralTransform(Nr,Nth,Nphi,Mr,Lmax,Lgrid,&
                                       &GLQWeights,GLQZeros,&
                                       &rho,theta,phi,&
                                       &dSdt,dadt)

        !$OMP PARALLEL DO
        DO n=1,(Mr+1)
            a(n,:,:,:) = 0.620101851488403D0 * aaRK(n,:,:,:,1) +&
                   & 0.379898148511597D0 * aaRK(n,:,:,:,3) +&
                   & 0.251891774271694D0 * dt * dadt(n,:,:,:)
        END DO
        !$OMP END PARALLEL DO


!--------------------------------------------------------!
!      STEP 4                                            !
!--------------------------------------------------------!

        !Store the coefficients inside aaRK 
        !which will be used in the Runge-Kutta method
        !$OMP PARALLEL DO
        DO n=1,(Mr+1)
            aaRK(n,:,:,:,4) = a(n,:,:,:)
        END DO
        !$OMP END PARALLEL DO

        CALL SpectralToSpatialTransform(Nr,Nth,Nphi,Mr,Lmax,Lgrid,&
                &GLQWeights,GLQZeros,rho,theta,phi,a,S)

        !Calculate the derivatives
        CALL EvaluatedSdr(Nr,Nth,Nphi,Mr,Lmax,Lgrid,GLQWeights,GLQZeros,&
                &rho,theta,phi,rmax,rmin,a,dSdr)
        CALL EvaluatedSdphi(Nr,Nth,Nphi,Mr,Lmax,Lgrid,GLQWeights,GLQZeros,&
                &rho,theta,phi,a,dSdphi)
        CALL EvaluatedSdtheta(Nr,Nth,Nphi,Mr,Lmax,Lgrid,GLQWeights,GLQZeros,&
                &rho,theta,phi,a,dSdth)
        CALL Evaluated2Sdr2(Nr,Nth,Nphi,Mr,Lmax,Lgrid,GLQWeights,GLQZeros,&
                &rho,theta,phi,rmax,rmin,dSdr,d2Sdr2)
        CALL Evaluated2Sdphi2(Nr,Nth,Nphi,Mr,Lmax,Lgrid,GLQWeights,GLQZeros,&
                &rho,theta,phi,dSdphi,d2Sdphi2)
        CALL Evaluated2Sdtheta2(Nr,Nth,Nphi,Mr,Lmax,Lgrid,GLQWeights,GLQZeros,&
                &rho,theta,phi,dSdth,d2Sdth2)
        
        !$OMP PARALLEL DO PRIVATE(j, k, rsqrd, bdryterm)
        DO i = 1,Nr
           DO j = 1,Nth
              DO k = 1,Nphi
                 
                 IF( i .EQ. 1 ) THEN
                    !BOUNDARIES (APPLY BOUNDARY CONDITIONS)
                    rsqrd = r(i) * r(i)
                    bdryterm = -1.0D0*epsStar*sigma*( S(i,j,k)**4 )

                    IF( theta(j) .LT. PI/2.0D0 ) THEN
                        bdryterm = bdryterm + &
                            &alpStar*SolPhi*COS(theta(j))
                    END IF

                    dSdt(i,j,k) = &
                        &2.0D0*bdryterm/r(i) - &
                        &4.0D0*epsStar*sigma*bdryterm*( S(i,j,k)**3 ) + &
                        &dSdth(i,j,k)/(rsqrd*TAN(theta(j))) + &
                        &d2Sdth2(i,j,k)/rsqrd + &
                        &d2Sdphi2(i,j,k)/(rsqrd*SIN(theta(j))*SIN(theta(j)))

                 ELSE
                    !INNER POINTS
                    rsqrd = r(i) * r(i)

                    dSdt(i,j,k) = &
                        &2.0D0*dSdr(i,j,k)/r(i) + d2Sdr2(i,j,k) + &
                        &dSdth(i,j,k)/(rsqrd*TAN(theta(j))) + &
                        &d2Sdth2(i,j,k)/rsqrd + &
                        &d2Sdphi2(i,j,k)/(rsqrd*SIN(theta(j))*SIN(theta(j)))

                 END IF

              END DO
           END DO
        END DO
        !$OMP END PARALLEL DO
 
        CALL SpatialToSpectralTransform(Nr,Nth,Nphi,Mr,Lmax,Lgrid,&
                                       &GLQWeights,GLQZeros,&
                                       &rho,theta,phi,&
                                       &dSdt,dadt)

        !$OMP PARALLEL DO
        DO n=1,(Mr+1)
            a(n,:,:,:) = 0.178079954393132D0 * aaRK(n,:,:,:,1) +& 
                   & 0.821920045606868D0 * aaRK(n,:,:,:,4) +&
                   & 0.544974750228521D0 * dt * dadt(n,:,:,:)
        END DO
        !$OMP END PARALLEL DO

!--------------------------------------------------------!
!      STEP 5                                            !
!--------------------------------------------------------!

        !Store the coefficients inside aaRK 
        !which will be used in the Runge-Kutta method
        !$OMP PARALLEL DO
        DO n=1,(Mr+1)
            aaRK(n,:,:,:,5) = a(n,:,:,:)
            dadt3(n,:,:,:) = dadt(n,:,:,:)
        END DO
        !$OMP END PARALLEL DO
 
        CALL SpectralToSpatialTransform(Nr,Nth,Nphi,Mr,Lmax,Lgrid,&
                &GLQWeights,GLQZeros,rho,theta,phi,a,S)

        !Calculate the derivatives
        CALL EvaluatedSdr(Nr,Nth,Nphi,Mr,Lmax,Lgrid,GLQWeights,GLQZeros,&
                &rho,theta,phi,rmax,rmin,a,dSdr)
        CALL EvaluatedSdphi(Nr,Nth,Nphi,Mr,Lmax,Lgrid,GLQWeights,GLQZeros,&
                &rho,theta,phi,a,dSdphi)
        CALL EvaluatedSdtheta(Nr,Nth,Nphi,Mr,Lmax,Lgrid,GLQWeights,GLQZeros,&
                &rho,theta,phi,a,dSdth)
        CALL Evaluated2Sdr2(Nr,Nth,Nphi,Mr,Lmax,Lgrid,GLQWeights,GLQZeros,&
                &rho,theta,phi,rmax,rmin,dSdr,d2Sdr2)
        CALL Evaluated2Sdphi2(Nr,Nth,Nphi,Mr,Lmax,Lgrid,GLQWeights,GLQZeros,&
                &rho,theta,phi,dSdphi,d2Sdphi2)
        CALL Evaluated2Sdtheta2(Nr,Nth,Nphi,Mr,Lmax,Lgrid,GLQWeights,GLQZeros,&
                &rho,theta,phi,dSdth,d2Sdth2)
        
        !$OMP PARALLEL DO PRIVATE(j, k, rsqrd, bdryterm)
        DO i = 1,Nr
           DO j = 1,Nth
              DO k = 1,Nphi
                 
                 IF( i .EQ. 1 ) THEN
                    !BOUNDARIES (APPLY BOUNDARY CONDITIONS)
                    rsqrd = r(i) * r(i)
                    bdryterm = -1.0D0*epsStar*sigma*( S(i,j,k)**4 )

                    IF( theta(j) .LT. PI/2.0D0 ) THEN
                        bdryterm = bdryterm + &
                            &alpStar*SolPhi*COS(theta(j))
                    END IF

                    dSdt(i,j,k) = &
                        &2.0D0*bdryterm/r(i) - &
                        &4.0D0*epsStar*sigma*bdryterm*( S(i,j,k)**3 ) + &
                        &dSdth(i,j,k)/(rsqrd*TAN(theta(j))) + &
                        &d2Sdth2(i,j,k)/rsqrd + &
                        &d2Sdphi2(i,j,k)/(rsqrd*SIN(theta(j))*SIN(theta(j)))

                 ELSE
                    !INNER POINTS
                    rsqrd = r(i) * r(i)

                    dSdt(i,j,k) = &
                        &2.0D0*dSdr(i,j,k)/r(i) + d2Sdr2(i,j,k) + &
                        &dSdth(i,j,k)/(rsqrd*TAN(theta(j))) + &
                        &d2Sdth2(i,j,k)/rsqrd + &
                        &d2Sdphi2(i,j,k)/(rsqrd*SIN(theta(j))*SIN(theta(j)))

                 END IF

              END DO
           END DO
        END DO
        !$OMP END PARALLEL DO

        CALL SpatialToSpectralTransform(Nr,Nth,Nphi,Mr,Lmax,Lgrid,&
                                       &GLQWeights,GLQZeros,&
                                       &rho,theta,phi,&
                                       &dSdt,dadt)

        !$OMP PARALLEL DO
        DO n=1,(Mr+1)
            a(n,:,:,:) = 0.517231671970585D0 * aaRK(n,:,:,:,3) + &
                   & 0.096059710526147D0 * aaRK(n,:,:,:,4) + &
                   & 0.063692468666290D0 * dt * dadt3(n,:,:,:) + &
                   & 0.386708617503269D0 * aaRK(n,:,:,:,5) + &
                   & 0.226007483236906D0 * dt * dadt(n,:,:,:)
        END DO
        !$OMP END PARALLEL DO

        PRINT *, 'time=', t
        PRINT *, 'dt=', dt
        t = t + dt
        
        RETURN
      END SUBROUTINE EvolveData
