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
&dt,&
&AF, AFinv, Dth, Dphi, F,&
&a)

        !EvolveData subroutine uses the Strong-Stability Preserving Runge-Kutta (SSPRK (3,3))
        !An outline of the integrator can be found on "Spectral Methods for Time-Dependent
        !Problems" by Jan Hesthaven, pg 201
        !u^(0) = a(t)
        !u^(1) = u^(0) + dt * L( u^(0) )
        !u^(2) = (3/4)u^(0) + (1/4)u^(1) + (1/4)dt * L( u^(1) )
        !a(t+dt) = u^(3) = (1/3)u^(0) + (2/3)u^(2) + (2/3)dt * L( u^(2) )
        !with CFL coefficient of 1
        
        USE             omp_lib
        IMPLICIT        none

!--------------------------------------------------------!
!     Declare calling variables                          !
!--------------------------------------------------------!

        INTEGER*4, INTENT(in)  :: M, Mr, NP, Lmax

        REAL*8, INTENT(in)   :: rootsign
        REAL*8, INTENT(in)   :: rmin
        REAL*8, INTENT(in)   :: rmax

        REAL*8, INTENT(in)   :: alpha(4*NP)
        REAL*8, INTENT(in)   :: betaR(4*NP)
        REAL*8, INTENT(in)   :: betaTh(4*NP)
        REAL*8, INTENT(in)   :: betaPhi(4*NP)

        REAL*8, INTENT(in)   :: gRR(4*NP)
        REAL*8, INTENT(in)   :: gThTh(4*NP)
        REAL*8, INTENT(in)   :: gPhiPhi(4*NP)

        REAL*8, INTENT(in)   :: gRTh(4*NP)
        REAL*8, INTENT(in)   :: gRPhi(4*NP)
        REAL*8, INTENT(in)   :: gThPhi(4*NP)

        REAL*8, INTENT(in)   :: dt
        REAL*8, INTENT(in)   :: F(NP)
        COMPLEX*16, INTENT(in)  :: AF(4*NP,NP)
        COMPLEX*16, INTENT(in)  :: AFinv(NP, 4*NP)
        COMPLEX*16, INTENT(in)  :: Dth(4*NP,NP)
        COMPLEX*16, INTENT(in)  :: Dphi(4*NP,NP)

        COMPLEX*16, INTENT(inout)  :: a(NP)

        COMPLEX*16      SphHarmonicY

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

        COMPLEX*16      aaRK(NP,3)
        COMPLEX*16      dadt(NP)

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

        a = aaRK(:,1) + dt*dadt


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

        a = 0.75D0*aaRK(:,1) + 0.25D0*aaRK(:,2) + 0.25D0*dt*dadt

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

        a = aaRK(:,1)/3.0D0 + 2.0D0*aaRK(:,3)/3.0D0 + 2.0D0*dt*dadt/3.0D0

        RETURN
      END SUBROUTINE EvolveData
