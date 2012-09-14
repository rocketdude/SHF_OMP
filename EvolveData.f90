!--------------------------------------------------------!
!    EvolveData subroutine                               !
!--------------------------------------------------------!

      SUBROUTINE EvolveData(&
&M, Mr, NP, Lmax,&
&rootsign, rmin, rmax,&
&alpha,&
&betaR, betaTh, betaPhi,&
&gRR, gThTh, gPhiPhi,&
&gRTh, gRPhi, gThPhi,&
&t, dt,&
&AF, AFinv, Dth, Dphi, F,&
&a)

        !EvolveData subroutine uses the Strong-Stability Preserving Runge-Kutta (SSPRK (5,4))
        !An outline of the integrator can be found on "Spectral Methods for Time-Dependent
        !Problems" by Jan Hesthaven, pg 201
        !with CFL coefficient of 1.508
        
        USE             omp_lib
        IMPLICIT        none

!--------------------------------------------------------!
!     Declare calling variables                          !
!--------------------------------------------------------!

        INTEGER*4               :: M, Mr, NP, Lmax

        REAL*8                  :: rootsign
        REAL*8                  :: rmin
        REAL*8                  :: rmax

        REAL*8                  :: alpha(4*NP)
        REAL*8                  :: betaR(4*NP)
        REAL*8                  :: betaTh(4*NP)
        REAL*8                  :: betaPhi(4*NP)

        REAL*8                  :: gRR(4*NP)
        REAL*8                  :: gThTh(4*NP)
        REAL*8                  :: gPhiPhi(4*NP)

        REAL*8                  :: gRTh(4*NP)
        REAL*8                  :: gRPhi(4*NP)
        REAL*8                  :: gThPhi(4*NP)

        REAL*8                  :: t
        REAL*8                  :: dt
        REAL*8                  :: F(NP)
        COMPLEX*16              :: AF(4*NP,NP)
        COMPLEX*16              :: AFinv(NP, 4*NP)
        COMPLEX*16              :: Dth(4*NP,NP)
        COMPLEX*16              :: Dphi(4*NP,NP)

        COMPLEX*16              :: a(NP)

        COMPLEX*16                 SphHarmonicY

!--------------------------------------------------------!
!     Declare Locals                                     !
!--------------------------------------------------------!

        CHARACTER*32    CTemp
        INTEGER*4         i, j, k
        INTEGER*4         ccol, crow
        
        COMPLEX*16      sqrtterm

        COMPLEX*16      dSdr(4*NP)
        COMPLEX*16      dSdth(4*NP)
        COMPLEX*16      dSdphi(4*NP)
        COMPLEX*16      dSdt(4*NP)

        COMPLEX*16      aaRK(NP,5)
        COMPLEX*16      dadt(NP)
        COMPLEX*16      dadt3(NP)

!--------------------------------------------------------!
!      Main Subroutine                                   !
!--------------------------------------------------------!

!--------------------------------------------------------!
!      STEP 1                                            !
!--------------------------------------------------------!

!!$        !Filter coefficients a with an exponential fitler
!!$        a = F*a

        !Store the coefficients inside aaRK which will be used in the Runge-Kutta method
        aaRK(:,1) = a

        !Calculate the derivatives S,r , S,theta, and S,phi
        CALL EvaluatedSdr(&
             &M, Mr, NP, Lmax,&
             &rmin, rmax,&
             &a, AF,&
             &dSdr)

        dSdphi = MATMUL(Dphi, a)
        dSdth = MATMUL(Dth, a)

        !$OMP PARALLEL DO PRIVATE(j, k, crow, sqrtterm)
        DO i = 1, (Mr+1)
           DO j = 1, 2*M
              DO k = 1, 2*M
                 
                 crow = (i-1)*(2*M)*(2*M) + (j-1)*(2*M) + k

                 IF( i .EQ. (Mr+1) ) THEN
                    !BOUNDARIES (APPLY BOUNDARY CONDITIONS)
                    sqrtterm = SQRT(&
                         &gThTh(crow)*dSdth(crow)*dSdth(crow) +&
                         &gPhiPhi(crow)*dSdphi(crow)*dSdphi(crow) +&
                         &2.0D0*gThPhi(crow)*dSdth(crow)*dSdphi(crow)&
                         &)

                    dSdt(crow) = &
                         &betaTh(crow)*dSdth(crow) + &
                         &betaPhi(crow)*dSdphi(crow) + &
                         &rootsign*alpha(crow)*sqrtterm
                 ELSE
                    !INNER POINTS
                    sqrtterm = SQRT(&
                         &gRR(crow)*dSdr(crow)*dSdr(crow) +&
                         &gThTh(crow)*dSdth(crow)*dSdth(crow) +&
                         &gPhiPhi(crow)*dSdphi(crow)*dSdphi(crow) +&
                         &2.0D0*gRTh(crow)*dSdr(crow)*dSdth(crow) +&
                         &2.0D0*gRPhi(crow)*dSdr(crow)*dSdphi(crow) +&
                         &2.0D0*gThPhi(crow)*dSdth(crow)*dSdphi(crow)&
                         &)

                    dSdt(crow) = &
                         &betaR(crow)*dSdr(crow) + &
                         &betaTh(crow)*dSdth(crow) + &
                         &betaPhi(crow)*dSdphi(crow) + &
                         &rootsign*alpha(crow)*sqrtterm
                 END IF

              END DO
           END DO
        END DO
        !$OMP END PARALLEL DO
        
        dadt = MATMUL( AFinv, dSdt)

        a = aaRK(:,1) + 0.391752226571890D0 * dt * dadt


!--------------------------------------------------------!
!      STEP 2                                            !
!--------------------------------------------------------!

!!$        !Filter coefficients a with an exponential fitler
!!$        a = F*a

        !Store the coefficients inside aaRK which will be used in the Runge-Kutta method
        aaRK(:,2) = a

        !Calculate the derivatives S,r , S,theta, and S,phi
        CALL EvaluatedSdr(&
             &M, Mr, NP, Lmax,&
             &rmin, rmax,&
             &a, AF,&
             &dSdr)

        dSdphi = MATMUL(Dphi, a)
        dSdth = MATMUL(Dth, a)

        !$OMP PARALLEL DO PRIVATE(j, k, crow, sqrtterm)
        DO i = 1, (Mr+1)
           DO j = 1, 2*M
              DO k = 1, 2*M
                 
                 crow = (i-1)*(2*M)*(2*M) + (j-1)*(2*M) + k

                 IF( i .EQ. (Mr+1) ) THEN
                    !BOUNDARIES (APPLY BOUNDARY CONDITIONS)
                    sqrtterm = SQRT(&
                         &gThTh(crow)*dSdth(crow)*dSdth(crow) +&
                         &gPhiPhi(crow)*dSdphi(crow)*dSdphi(crow) +&
                         &2.0D0*gThPhi(crow)*dSdth(crow)*dSdphi(crow)&
                         &)

                    dSdt(crow) = &
                         &betaTh(crow)*dSdth(crow) + &
                         &betaPhi(crow)*dSdphi(crow) + &
                         &rootsign*alpha(crow)*sqrtterm
                 ELSE
                    !INNER POINTS
                    sqrtterm = SQRT(&
                         &gRR(crow)*dSdr(crow)*dSdr(crow) +&
                         &gThTh(crow)*dSdth(crow)*dSdth(crow) +&
                         &gPhiPhi(crow)*dSdphi(crow)*dSdphi(crow) +&
                         &2.0D0*gRTh(crow)*dSdr(crow)*dSdth(crow) +&
                         &2.0D0*gRPhi(crow)*dSdr(crow)*dSdphi(crow) +&
                         &2.0D0*gThPhi(crow)*dSdth(crow)*dSdphi(crow)&
                         &)

                    dSdt(crow) = &
                         &betaR(crow)*dSdr(crow) + &
                         &betaTh(crow)*dSdth(crow) + &
                         &betaPhi(crow)*dSdphi(crow) + &
                         &rootsign*alpha(crow)*sqrtterm
                 END IF

              END DO
           END DO
        END DO
        !$OMP END PARALLEL DO
        
        dadt = MATMUL( AFinv, dSdt)

        a = 0.444370493651235D0 * aaRK(:,1) + 0.555629506348765D0 * aaRK(:,2)&
            &+ 0.368410593050371D0 * dt * dadt

!--------------------------------------------------------!
!      STEP 3                                            !
!--------------------------------------------------------!

!!$        !Filter coefficients a with an exponential fitler
!!$        a = F*a

        !Store the coefficients inside aaRK which will be used in the Runge-Kutta method
        aaRK(:,3) = a

        !Calculate the derivatives S,r , S,theta, and S,phi
        CALL EvaluatedSdr(&
             &M, Mr, NP, Lmax,&
             &rmin, rmax,&
             &a, AF,&
             &dSdr)

        dSdphi = MATMUL(Dphi, a)
        dSdth = MATMUL(Dth, a)

        !$OMP PARALLEL DO PRIVATE(j, k, crow, sqrtterm)
        DO i = 1, (Mr+1)
           DO j = 1, 2*M
              DO k = 1, 2*M
                 
                 crow = (i-1)*(2*M)*(2*M) + (j-1)*(2*M) + k

                 IF( i .EQ. (Mr+1) ) THEN
                    !BOUNDARIES (APPLY BOUNDARY CONDITIONS)
                    sqrtterm = SQRT(&
                         &gThTh(crow)*dSdth(crow)*dSdth(crow) +&
                         &gPhiPhi(crow)*dSdphi(crow)*dSdphi(crow) +&
                         &2.0D0*gThPhi(crow)*dSdth(crow)*dSdphi(crow)&
                         &)

                    dSdt(crow) = &
                         &betaTh(crow)*dSdth(crow) + &
                         &betaPhi(crow)*dSdphi(crow) + &
                         &rootsign*alpha(crow)*sqrtterm
                 ELSE
                    !INNER POINTS
                    sqrtterm = SQRT(&
                         &gRR(crow)*dSdr(crow)*dSdr(crow) +&
                         &gThTh(crow)*dSdth(crow)*dSdth(crow) +&
                         &gPhiPhi(crow)*dSdphi(crow)*dSdphi(crow) +&
                         &2.0D0*gRTh(crow)*dSdr(crow)*dSdth(crow) +&
                         &2.0D0*gRPhi(crow)*dSdr(crow)*dSdphi(crow) +&
                         &2.0D0*gThPhi(crow)*dSdth(crow)*dSdphi(crow)&
                         &)

                    dSdt(crow) = &
                         &betaR(crow)*dSdr(crow) + &
                         &betaTh(crow)*dSdth(crow) + &
                         &betaPhi(crow)*dSdphi(crow) + &
                         &rootsign*alpha(crow)*sqrtterm
                 END IF

              END DO
           END DO
        END DO
        !$OMP END PARALLEL DO
        
        dadt = MATMUL( AFinv, dSdt)

        a = 0.620101851488403D0 * aaRK(:,1) + 0.379898148511597D0 * aaRK(:,3)&
            & + 0.251891774271694D0 * dt * dadt

!--------------------------------------------------------!
!      STEP 4                                            !
!--------------------------------------------------------!

!!$        !Filter coefficients a with an exponential fitler
!!$        a = F*a

        !Store the coefficients inside aaRK which will be used in the Runge-Kutta method
        aaRK(:,4) = a

        !Calculate the derivatives S,r , S,theta, and S,phi
        CALL EvaluatedSdr(&
             &M, Mr, NP, Lmax,&
             &rmin, rmax,&
             &a, AF,&
             &dSdr)

        dSdphi = MATMUL(Dphi, a)
        dSdth = MATMUL(Dth, a)

        !$OMP PARALLEL DO PRIVATE(j, k, crow, sqrtterm)
        DO i = 1, (Mr+1)
           DO j = 1, 2*M
              DO k = 1, 2*M
                 
                 crow = (i-1)*(2*M)*(2*M) + (j-1)*(2*M) + k

                 IF( i .EQ. (Mr+1) ) THEN
                    !BOUNDARIES (APPLY BOUNDARY CONDITIONS)
                    sqrtterm = SQRT(&
                         &gThTh(crow)*dSdth(crow)*dSdth(crow) +&
                         &gPhiPhi(crow)*dSdphi(crow)*dSdphi(crow) +&
                         &2.0D0*gThPhi(crow)*dSdth(crow)*dSdphi(crow)&
                         &)

                    dSdt(crow) = &
                         &betaTh(crow)*dSdth(crow) + &
                         &betaPhi(crow)*dSdphi(crow) + &
                         &rootsign*alpha(crow)*sqrtterm
                 ELSE
                    !INNER POINTS
                    sqrtterm = SQRT(&
                         &gRR(crow)*dSdr(crow)*dSdr(crow) +&
                         &gThTh(crow)*dSdth(crow)*dSdth(crow) +&
                         &gPhiPhi(crow)*dSdphi(crow)*dSdphi(crow) +&
                         &2.0D0*gRTh(crow)*dSdr(crow)*dSdth(crow) +&
                         &2.0D0*gRPhi(crow)*dSdr(crow)*dSdphi(crow) +&
                         &2.0D0*gThPhi(crow)*dSdth(crow)*dSdphi(crow)&
                         &)

                    dSdt(crow) = &
                         &betaR(crow)*dSdr(crow) + &
                         &betaTh(crow)*dSdth(crow) + &
                         &betaPhi(crow)*dSdphi(crow) + &
                         &rootsign*alpha(crow)*sqrtterm
                 END IF

              END DO
           END DO
        END DO
        !$OMP END PARALLEL DO
        
        dadt = MATMUL( AFinv, dSdt)

        a = 0.178079954393132D0 * aaRK(:,1) + 0.821920045606868D0 * aaRK(:,4)&
            & + 0.544974750228521D0 * dt * dadt


!--------------------------------------------------------!
!      STEP 5                                            !
!--------------------------------------------------------!

!!$        !Filter coefficients a with an exponential fitler
!!$        a = F*a

        !Store the coefficients inside aaRK which will be used in the Runge-Kutta method
        aaRK(:,5) = a
        dadt3 = dadt

        !Calculate the derivatives S,r , S,theta, and S,phi
        CALL EvaluatedSdr(&
             &M, Mr, NP, Lmax,&
             &rmin, rmax,&
             &a, AF,&
             &dSdr)

        dSdphi = MATMUL(Dphi, a)
        dSdth = MATMUL(Dth, a)

        !$OMP PARALLEL DO PRIVATE(j, k, crow, sqrtterm)
        DO i = 1, (Mr+1)
           DO j = 1, 2*M
              DO k = 1, 2*M
                 
                 crow = (i-1)*(2*M)*(2*M) + (j-1)*(2*M) + k

                 IF( i .EQ. (Mr+1) ) THEN
                    !BOUNDARIES (APPLY BOUNDARY CONDITIONS)
                    sqrtterm = SQRT(&
                         &gThTh(crow)*dSdth(crow)*dSdth(crow) +&
                         &gPhiPhi(crow)*dSdphi(crow)*dSdphi(crow) +&
                         &2.0D0*gThPhi(crow)*dSdth(crow)*dSdphi(crow)&
                         &)

                    dSdt(crow) = &
                         &betaTh(crow)*dSdth(crow) + &
                         &betaPhi(crow)*dSdphi(crow) + &
                         &rootsign*alpha(crow)*sqrtterm
                 ELSE
                    !INNER POINTS
                    sqrtterm = SQRT(&
                         &gRR(crow)*dSdr(crow)*dSdr(crow) +&
                         &gThTh(crow)*dSdth(crow)*dSdth(crow) +&
                         &gPhiPhi(crow)*dSdphi(crow)*dSdphi(crow) +&
                         &2.0D0*gRTh(crow)*dSdr(crow)*dSdth(crow) +&
                         &2.0D0*gRPhi(crow)*dSdr(crow)*dSdphi(crow) +&
                         &2.0D0*gThPhi(crow)*dSdth(crow)*dSdphi(crow)&
                         &)

                    dSdt(crow) = &
                         &betaR(crow)*dSdr(crow) + &
                         &betaTh(crow)*dSdth(crow) + &
                         &betaPhi(crow)*dSdphi(crow) + &
                         &rootsign*alpha(crow)*sqrtterm
                 END IF

              END DO
           END DO
        END DO
        !$OMP END PARALLEL DO
        
        dadt = MATMUL( AFinv, dSdt)

        a = 0.517231671970585D0 * aaRK(:,3) + 0.096059710526147D0 * aaRK(:,4)&
            & + 0.063692468666290D0 * dt * dadt3 + 0.386708617503269D0 * aaRK(:,5)&
            & + 0.226007483236906D0 * dt * dadt

        t = t + dt
        
        RETURN
      END SUBROUTINE EvolveData
