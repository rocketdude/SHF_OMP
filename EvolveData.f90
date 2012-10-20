!--------------------------------------------------------!
!    EvolveData subroutine                               !
!--------------------------------------------------------!

      SUBROUTINE EvolveData(&
&Nr, Nth, Nphi, Mr, Lmax, Lgrid,&
&GLQWeights, GLQZeros,&
&rootsign, rmin, rmax,&
&rho, theta, phi,&
&alpha,&
&betaR, betaTh, betaPhi,&
&gRR, gThTh, gPhiPhi,&
&gRTh, gRPhi, gThPhi,&
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

        REAL*8                  :: rootsign
        REAL*8                  :: rmin
        REAL*8                  :: rmax
        REAL*8                  :: GLQWeights(Lgrid+1), GLQZeros(Lgrid+1)
        REAL*8                  :: rho(Nr), theta(Nth), phi(Nphi)

        REAL*8                  :: alpha(Nr,Nth,Nphi)
        REAL*8                  :: betaR(Nr,Nth,Nphi)
        REAL*8                  :: betaTh(Nr,Nth,Nphi)
        REAL*8                  :: betaPhi(Nr,Nth,Nphi)

        REAL*8                  :: gRR(Nr,Nth,Nphi)
        REAL*8                  :: gThTh(Nr,Nth,Nphi)
        REAL*8                  :: gPhiPhi(Nr,Nth,Nphi)

        REAL*8                  :: gRTh(Nr,Nth,Nphi)
        REAL*8                  :: gRPhi(Nr,Nth,Nphi)
        REAL*8                  :: gThPhi(Nr,Nth,Nphi)

        REAL*8                  :: t
        REAL*8                  :: dt

        COMPLEX*16              :: a(Mr+1,2,Lmax+1,Lmax+1)

!--------------------------------------------------------!
!     Declare Locals                                     !
!--------------------------------------------------------!

        CHARACTER*32    CTemp
        INTEGER*4         i, j, k
        INTEGER*4         n
        
        COMPLEX*16      sqrtterm

        COMPLEX*16      dSdr(Nr,Nth,Nphi)
        COMPLEX*16      dSdth(Nr,Nth,Nphi)
        COMPLEX*16      dSdphi(Nr,Nth,Nphi)
        COMPLEX*16      dSdt(Nr,Nth,Nphi)

        COMPLEX*16      aaRK(Mr+1,2,Lmax+1,Lmax+1,5)
        COMPLEX*16      dadt(Mr+1,2,Lmax+1,Lmax+1)
        COMPLEX*16      dadt3(Mr+1,2,Lmax+1,Lmax+1)

!--------------------------------------------------------!
!      Main Subroutine                                   !
!--------------------------------------------------------!

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

        !Calculate the derivatives S,r , S,theta, and S,phi
        CALL EvaluatedSdr(Nr,Nth,Nphi,Mr,Lmax,Lgrid,GLQWeights,GLQZeros,&
                &rho,theta,phi,rmax,rmin,a,dSdr)
        CALL EvaluatedSdphi(Nr,Nth,Nphi,Mr,Lmax,Lgrid,GLQWeights,GLQZeros,&
                &rho,theta,phi,a,dSdphi)
        CALL EvaluatedSdtheta(Nr,Nth,Nphi,Mr,Lmax,Lgrid,GLQWeights,GLQZeros,&
                &rho,theta,phi,a,dSdth)
        
        !$OMP PARALLEL DO PRIVATE(j, k, sqrtterm)
        DO i = 1,Nr
           DO j = 1,Nth
              DO k = 1,Nphi
                 
                 IF( i .EQ. (Mr+1) ) THEN
                    !BOUNDARIES (APPLY BOUNDARY CONDITIONS)
                    sqrtterm = SQRT(&
                         &gThTh(i,j,k)*dSdth(i,j,k)*dSdth(i,j,k) +&
                         &gPhiPhi(i,j,k)*dSdphi(i,j,k)*dSdphi(i,j,k) +&
                         &2.0D0*gThPhi(i,j,k)*dSdth(i,j,k)*dSdphi(i,j,k)&
                         &)

                    dSdt(i,j,k) = &
                         &betaTh(i,j,k)*dSdth(i,j,k) + &
                         &betaPhi(i,j,k)*dSdphi(i,j,k) + &
                         &rootsign*alpha(i,j,k)*sqrtterm
                 ELSE
                    !INNER POINTS
                    sqrtterm = SQRT(&
                         &gRR(i,j,k)*dSdr(i,j,k)*dSdr(i,j,k) +&
                         &gThTh(i,j,k)*dSdth(i,j,k)*dSdth(i,j,k) +&
                         &gPhiPhi(i,j,k)*dSdphi(i,j,k)*dSdphi(i,j,k) +&
                         &2.0D0*gRTh(i,j,k)*dSdr(i,j,k)*dSdth(i,j,k) +&
                         &2.0D0*gRPhi(i,j,k)*dSdr(i,j,k)*dSdphi(i,j,k) +&
                         &2.0D0*gThPhi(i,j,k)*dSdth(i,j,k)*dSdphi(i,j,k)&
                         &)

                    dSdt(i,j,k) = &
                         &betaR(i,j,k)*dSdr(i,j,k) + &
                         &betaTh(i,j,k)*dSdth(i,j,k) + &
                         &betaPhi(i,j,k)*dSdphi(i,j,k) + &
                         &rootsign*alpha(i,j,k)*sqrtterm
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

        !Calculate the derivatives S,r , S,theta, and S,phi
        CALL EvaluatedSdr(Nr,Nth,Nphi,Mr,Lmax,Lgrid,GLQWeights,GLQZeros,&
                &rho,theta,phi,rmax,rmin,a,dSdr)
        CALL EvaluatedSdphi(Nr,Nth,Nphi,Mr,Lmax,Lgrid,GLQWeights,GLQZeros,&
                &rho,theta,phi,a,dSdphi)
        CALL EvaluatedSdtheta(Nr,Nth,Nphi,Mr,Lmax,Lgrid,GLQWeights,GLQZeros,&
                &rho,theta,phi,a,dSdth)

        !$OMP PARALLEL DO PRIVATE(j, k, sqrtterm)
        DO i = 1,Nr
           DO j = 1,Nth
              DO k = 1,Nphi
                 
                 IF( i .EQ. (Mr+1) ) THEN
                    !BOUNDARIES (APPLY BOUNDARY CONDITIONS)
                    sqrtterm = SQRT(&
                         &gThTh(i,j,k)*dSdth(i,j,k)*dSdth(i,j,k) +&
                         &gPhiPhi(i,j,k)*dSdphi(i,j,k)*dSdphi(i,j,k) +&
                         &2.0D0*gThPhi(i,j,k)*dSdth(i,j,k)*dSdphi(i,j,k)&
                         &)

                    dSdt(i,j,k) = &
                         &betaTh(i,j,k)*dSdth(i,j,k) + &
                         &betaPhi(i,j,k)*dSdphi(i,j,k) + &
                         &rootsign*alpha(i,j,k)*sqrtterm
                 ELSE
                    !INNER POINTS
                    sqrtterm = SQRT(&
                         &gRR(i,j,k)*dSdr(i,j,k)*dSdr(i,j,k) +&
                         &gThTh(i,j,k)*dSdth(i,j,k)*dSdth(i,j,k) +&
                         &gPhiPhi(i,j,k)*dSdphi(i,j,k)*dSdphi(i,j,k) +&
                         &2.0D0*gRTh(i,j,k)*dSdr(i,j,k)*dSdth(i,j,k) +&
                         &2.0D0*gRPhi(i,j,k)*dSdr(i,j,k)*dSdphi(i,j,k) +&
                         &2.0D0*gThPhi(i,j,k)*dSdth(i,j,k)*dSdphi(i,j,k)&
                         &)

                    dSdt(i,j,k) = &
                         &betaR(i,j,k)*dSdr(i,j,k) + &
                         &betaTh(i,j,k)*dSdth(i,j,k) + &
                         &betaPhi(i,j,k)*dSdphi(i,j,k) + &
                         &rootsign*alpha(i,j,k)*sqrtterm
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

        !Calculate the derivatives S,r , S,theta, and S,phi
        CALL EvaluatedSdr(Nr,Nth,Nphi,Mr,Lmax,Lgrid,GLQWeights,GLQZeros,&
                &rho,theta,phi,rmax,rmin,a,dSdr)
        CALL EvaluatedSdphi(Nr,Nth,Nphi,Mr,Lmax,Lgrid,GLQWeights,GLQZeros,&
                &rho,theta,phi,a,dSdphi)
        CALL EvaluatedSdtheta(Nr,Nth,Nphi,Mr,Lmax,Lgrid,GLQWeights,GLQZeros,&
                &rho,theta,phi,a,dSdth)

        !$OMP PARALLEL DO PRIVATE(j, k, sqrtterm)
        DO i = 1,Nr
           DO j = 1,Nth
              DO k = 1,Nphi
                 
                 IF( i .EQ. (Mr+1) ) THEN
                    !BOUNDARIES (APPLY BOUNDARY CONDITIONS)
                    sqrtterm = SQRT(&
                         &gThTh(i,j,k)*dSdth(i,j,k)*dSdth(i,j,k) +&
                         &gPhiPhi(i,j,k)*dSdphi(i,j,k)*dSdphi(i,j,k) +&
                         &2.0D0*gThPhi(i,j,k)*dSdth(i,j,k)*dSdphi(i,j,k)&
                         &)

                    dSdt(i,j,k) = &
                         &betaTh(i,j,k)*dSdth(i,j,k) + &
                         &betaPhi(i,j,k)*dSdphi(i,j,k) + &
                         &rootsign*alpha(i,j,k)*sqrtterm
                 ELSE
                    !INNER POINTS
                    sqrtterm = SQRT(&
                         &gRR(i,j,k)*dSdr(i,j,k)*dSdr(i,j,k) +&
                         &gThTh(i,j,k)*dSdth(i,j,k)*dSdth(i,j,k) +&
                         &gPhiPhi(i,j,k)*dSdphi(i,j,k)*dSdphi(i,j,k) +&
                         &2.0D0*gRTh(i,j,k)*dSdr(i,j,k)*dSdth(i,j,k) +&
                         &2.0D0*gRPhi(i,j,k)*dSdr(i,j,k)*dSdphi(i,j,k) +&
                         &2.0D0*gThPhi(i,j,k)*dSdth(i,j,k)*dSdphi(i,j,k)&
                         &)

                    dSdt(i,j,k) = &
                         &betaR(i,j,k)*dSdr(i,j,k) + &
                         &betaTh(i,j,k)*dSdth(i,j,k) + &
                         &betaPhi(i,j,k)*dSdphi(i,j,k) + &
                         &rootsign*alpha(i,j,k)*sqrtterm
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

        !Calculate the derivatives S,r , S,theta, and S,phi
        CALL EvaluatedSdr(Nr,Nth,Nphi,Mr,Lmax,Lgrid,GLQWeights,GLQZeros,&
                &rho,theta,phi,rmax,rmin,a,dSdr)
        CALL EvaluatedSdphi(Nr,Nth,Nphi,Mr,Lmax,Lgrid,GLQWeights,GLQZeros,&
                &rho,theta,phi,a,dSdphi)
        CALL EvaluatedSdtheta(Nr,Nth,Nphi,Mr,Lmax,Lgrid,GLQWeights,GLQZeros,&
                &rho,theta,phi,a,dSdth)

        !$OMP PARALLEL DO PRIVATE(j, k, sqrtterm)
        DO i = 1,Nr
           DO j = 1,Nth
              DO k = 1,Nphi
                 
                 IF( i .EQ. (Mr+1) ) THEN
                    !BOUNDARIES (APPLY BOUNDARY CONDITIONS)
                    sqrtterm = SQRT(&
                         &gThTh(i,j,k)*dSdth(i,j,k)*dSdth(i,j,k) +&
                         &gPhiPhi(i,j,k)*dSdphi(i,j,k)*dSdphi(i,j,k) +&
                         &2.0D0*gThPhi(i,j,k)*dSdth(i,j,k)*dSdphi(i,j,k)&
                         &)

                    dSdt(i,j,k) = &
                         &betaTh(i,j,k)*dSdth(i,j,k) + &
                         &betaPhi(i,j,k)*dSdphi(i,j,k) + &
                         &rootsign*alpha(i,j,k)*sqrtterm
                 ELSE
                    !INNER POINTS
                    sqrtterm = SQRT(&
                         &gRR(i,j,k)*dSdr(i,j,k)*dSdr(i,j,k) +&
                         &gThTh(i,j,k)*dSdth(i,j,k)*dSdth(i,j,k) +&
                         &gPhiPhi(i,j,k)*dSdphi(i,j,k)*dSdphi(i,j,k) +&
                         &2.0D0*gRTh(i,j,k)*dSdr(i,j,k)*dSdth(i,j,k) +&
                         &2.0D0*gRPhi(i,j,k)*dSdr(i,j,k)*dSdphi(i,j,k) +&
                         &2.0D0*gThPhi(i,j,k)*dSdth(i,j,k)*dSdphi(i,j,k)&
                         &)

                    dSdt(i,j,k) = &
                         &betaR(i,j,k)*dSdr(i,j,k) + &
                         &betaTh(i,j,k)*dSdth(i,j,k) + &
                         &betaPhi(i,j,k)*dSdphi(i,j,k) + &
                         &rootsign*alpha(i,j,k)*sqrtterm
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

        !Calculate the derivatives S,r , S,theta, and S,phi
        CALL EvaluatedSdr(Nr,Nth,Nphi,Mr,Lmax,Lgrid,GLQWeights,GLQZeros,&
                &rho,theta,phi,rmax,rmin,a,dSdr)
        CALL EvaluatedSdphi(Nr,Nth,Nphi,Mr,Lmax,Lgrid,GLQWeights,GLQZeros,&
                &rho,theta,phi,a,dSdphi)
        CALL EvaluatedSdtheta(Nr,Nth,Nphi,Mr,Lmax,Lgrid,GLQWeights,GLQZeros,&
                &rho,theta,phi,a,dSdth)

        !$OMP PARALLEL DO PRIVATE(j, k, sqrtterm)
        DO i = 1,Nr
           DO j = 1,Nth
              DO k = 1,Nphi
                 
                 IF( i .EQ. (Mr+1) ) THEN
                    !BOUNDARIES (APPLY BOUNDARY CONDITIONS)
                    sqrtterm = SQRT(&
                         &gThTh(i,j,k)*dSdth(i,j,k)*dSdth(i,j,k) +&
                         &gPhiPhi(i,j,k)*dSdphi(i,j,k)*dSdphi(i,j,k) +&
                         &2.0D0*gThPhi(i,j,k)*dSdth(i,j,k)*dSdphi(i,j,k)&
                         &)

                    dSdt(i,j,k) = &
                         &betaTh(i,j,k)*dSdth(i,j,k) + &
                         &betaPhi(i,j,k)*dSdphi(i,j,k) + &
                         &rootsign*alpha(i,j,k)*sqrtterm
                 ELSE
                    !INNER POINTS
                    sqrtterm = SQRT(&
                         &gRR(i,j,k)*dSdr(i,j,k)*dSdr(i,j,k) +&
                         &gThTh(i,j,k)*dSdth(i,j,k)*dSdth(i,j,k) +&
                         &gPhiPhi(i,j,k)*dSdphi(i,j,k)*dSdphi(i,j,k) +&
                         &2.0D0*gRTh(i,j,k)*dSdr(i,j,k)*dSdth(i,j,k) +&
                         &2.0D0*gRPhi(i,j,k)*dSdr(i,j,k)*dSdphi(i,j,k) +&
                         &2.0D0*gThPhi(i,j,k)*dSdth(i,j,k)*dSdphi(i,j,k)&
                         &)

                    dSdt(i,j,k) = &
                         &betaR(i,j,k)*dSdr(i,j,k) + &
                         &betaTh(i,j,k)*dSdth(i,j,k) + &
                         &betaPhi(i,j,k)*dSdphi(i,j,k) + &
                         &rootsign*alpha(i,j,k)*sqrtterm
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

        t = t + dt
        
        RETURN
      END SUBROUTINE EvolveData
