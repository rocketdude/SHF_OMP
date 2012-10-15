!========================================================!
!    FindU Subroutine                                    !
!========================================================!

     SUBROUTINE FindU(&
&Nr, Nth, Nphi, Mr, Mlm, SpM,&
&gRR, gThTh, gPhiPhi,&
&gRTh, gRPhi, gThPhi,&
&r, rho, theta, phi,&
&a,&
&it, WriteSit,&
&U, Uave, USp,&
&thetaSp, phiSp,&
&gRRUsqrd, gRRUsqrdAve)

!This subroutine finds U(theta, phi) using cubic spline interpolation
!and then performing linear interpolation

        USE             omp_lib
        IMPLICIT        none

!-----------------------------------------------------------!
!       Declare calling variables                           !
!-----------------------------------------------------------!

       INTEGER*4            :: Nr, Nth, Nphi, Mr, Mlm, SpM, it, WriteSit

       REAL*8               :: r(Nr), rho(Nr), theta(Nth), phi(Nphi)

       REAL*8               :: thetaSp(SpM)
       REAL*8               :: phiSp(SpM)

       REAL*8               :: gRR(Nr,Nth,Nphi)
       REAL*8               :: gThTh(Nr,Nth,Nphi)
       REAL*8               :: gPhiPhi(Nr,Nth,Nphi)
       REAL*8               :: gRTh(Nr,Nth,Nphi)
       REAL*8               :: gRPhi(Nr,Nth,Nphi)
       REAL*8               :: gThPhi(Nr,Nth,Nphi)
    
       COMPLEX*16           :: a(Mr+1, Mlm)

       REAL*8, INTENT(out)  :: U(Nth, Nphi)
       REAL*8, INTENT(out)  :: Uave
       REAL*8, INTENT(out)  :: USp(SpM)
       REAL*8, INTENT(out)  :: gRRUsqrd(Nth, Nphi)
       REAL*8, INTENT(out)  :: gRRUsqrdAve

!--------------------------------------------------------!
!     Declare Locals                                     !
!--------------------------------------------------------!

       CHARACTER*32      CTemp

       INTEGER*4, PARAMETER ::  SP=1000 
       !Number of points to be calculated on the spline
       INTEGER*4, PARAMETER ::  PolyInterpOrder=6 
       !Number of points for poly. interp.
       INTEGER*4                i, j, k, ii, dii, jj, ss
       INTEGER*4                n, l, ml
       INTEGER*4                crow
       INTEGER*4                iindex, jkindex
       INTEGER*4                ilow, ilow2, flag
       INTEGER*4                error

       COMPLEX*16             S(Nr,Nth,Nphi)
       COMPLEX*16             TnYlm
    
       REAL*8                 U_r(Nr) 
       REAL*8                 U_r2(Nr) 
                              !Second derivatives of U_r at points along r

       REAL*8                 g_r(6,Nr)
       REAL*8                 gUU(6)
       REAL*8                 gDD(6)
       REAL*8                 gPolyInt(6,PolyInterpOrder)
       REAL*8                 rPolyInt(PolyInterpOrder)

       REAL*8                 rr(SP)
       REAL*8                 UU(SP)
       REAL*8                 Temp
       
       REAL*8                 TargetedSValue
       REAL*8                 deltar

!--------------------------------------------------------!
!      Main Subroutine                                   !
!--------------------------------------------------------!

       TargetedSValue = 100.0D0

       !Evaluate the eikonal data S
       CALL SpectralToSpatialTransform(Nr, Nth, Nphi, Mr, Mlm,&
                                      &rho, theta, phi, a, S)

       !Writing S into file (only at certain iterations)
       IF( MOD(it, WriteSit) .EQ. 0 ) CALL WriteS(Nr, Nth, Nphi, ABS(S), it)

       !$OMP PARALLEL DO &
       !$OMP &PRIVATE(k, i, iindex, jkindex, U_r, U_r2, flag, gDD, gUU,&
       !$OMP &       g_r, rPolyInt, gPolyInt, jj, dii, ilow, ilow2,&
       !$OMP &       deltar, rr, UU)
       DO j = 1, Nth
          DO k = 1, Nphi

             U_r = ABS(S(:,j,k))

             !To find the radial distance of the light cone,
             !we fit cubic splines and then perform linear interpolation
             
             !Calculate the second derivatives and store it into U_r2

             CALL ComputeSpline2ndDeriv(r, U_r, Nr, 1.0D31, 1.0D31, U_r2)
                    !Using natural cubic spline

             !Find ilow and ilow+1 which are 
             !the indices that bound where TargetedSValue is located
             flag = 0
             DO i=(1+(Nr/10)), (Nr-(Nr/10))
                IF( U_r(i) .GE. TargetedSValue ) THEN
                   ilow = i-1
                   flag = 1
                   EXIT
                END IF
             END DO

             IF( flag .EQ. 0) THEN           
                WRITE(*,*) 'Error in hunting'
                !READ(*, '()')
                STOP
             END IF

             !Calculate the equispaced points rr--between r(ilow) and r(ilow+1)
             !and, calculate the values of ABS(S) at those points
             deltar = (r(ilow+1) - r(ilow))/DBLE(SP-1)
             DO i = 1, SP
                rr(i) = r(ilow) + DBLE(i-1)*deltar
                CALL CubicSplineInterpolation(r, U_r, U_r2, ilow, Nr,&
                                             &rr(i), UU(i))

             END DO

             !Perform HUNT again,finding ilow2 and ilow2+1 which are 
             !the indices that bound where TargetedSValue is located
             flag = 0
             DO i=1, SP
                IF( UU(i) .GE. TargetedSValue ) THEN
                   ilow2 = i-1
                   flag = 1
                   EXIT
                END IF
             END DO

             IF( flag .EQ. 0) THEN           
                WRITE(*,*) 'Error in hunting'
                !READ(*, '()')
                STOP
             END IF

             !Find the value of r where TargetedSValue is located by interp.
             U(j,k) = rr(ilow2) + &
                  &( rr(ilow2+1) - rr(ilow2) )/( UU(ilow2+1) - UU(ilow2) ) *&
                  &(TargetedSValue - UU(ilow2) )
    
             g_r(1,:) = gRR(:,j,k)
             g_r(2,:) = gThTh(:,j,k)
             g_r(3,:) = gPhiPhi(:,j,k)
             g_r(4,:) = gRTh(:,j,k)
             g_r(5,:) = gRPhi(:,j,k)
             g_r(6,:) = gThPhi(:,j,k)
            
             !Interpolate to find the three metric at U(j,k)
             dii = PolyInterpOrder/2 !This is floored
             DO ii = 1, PolyInterpOrder
                rPolyInt(ii) = r(ilow - dii + ii)

                DO jj = 1, 6 !there are only 6 independent 3-metric components
                    gPolyInt(jj,ii) = g_r(jj,ilow - dii + ii)
                END DO
                
             END DO
             
             DO jj = 1,6
                CALL PolynomialInterpolation(rPolyInt, gPolyInt(jj,:),&
                                            &PolyInterpOrder, U(j,k), gUU(jj))
             END DO

             !Invert to get the lower index 3-metric
             CALL Invert3Metric(gUU(1), gUU(2), gUU(3),&
                               &gUU(4), gUU(5), gUU(6),&
                               &gDD(1), gDD(2), gDD(3),&
                               &gDD(4), gDD(5), gDD(6),&
                               &error) 
             
             !Now calculate g_rr * U^2
             gRRUsqrd(j,k) = gDD(1) * U(j,k) * U(j,k)
             
          END DO
       END DO
       !$OMP END PARALLEL DO

       !Calculate the average U: Uave 
       !(only applies if there's obvious spherical symmetry)
       Uave = SUM( U )/ (Nth*Nphi)
       !Calculate the average g_rr*U^2 (the above premise applies here as well)
       gRRUsqrdAve = SUM( gRRUsqrd ) / (Nth*Nphi)

       !Calculate the values of U at the requested directions
       !$OMP PARALLEL DO &
       !$OMP &PRIVATE(i, n, l, ml, U_r, U_r2, flag,&
       !$OMP &        ilow, ilow2, deltar, rr, UU, crow, TnYlm)
       DO ss = 1, SpM

             CALL SpectralToSpatialOnARadialLine(Nr,Nth,Nphi,Mr,Mlm,&
                                                &rho,theta,phi,&
                                                &thetaSp(ss), phiSp(ss),&
                                                &a, U_r)
                            
             !Calculate the second derivatives and store it into U_r2
             CALL ComputeSpline2ndDeriv(r, U_r, Nr, 1.0D31, 1.0D31, U_r2)

             !Find ilow and ilow+1 which are the indices that 
             !bound where TargetedSValue is located
             flag = 0
             DO i=(1+(Nr/10)), (Nr-(Nr/10))
                IF( U_r(i) .GE. TargetedSValue ) THEN
                   ilow = i-1
                   flag = 1
                   EXIT
                END IF
             END DO

             IF( flag .EQ. 0) THEN           
                WRITE(*,*) 'Error in hunting for special points'
                !READ(*, '()')
                STOP
             END IF

             !Calculate the equispaced points rr--between r(ilow) and r(ilow+1)
             !and, calculate the values of ABS(S) at those points
             deltar = (r(ilow+1) - r(ilow))/DBLE(SP-1)
             DO i = 1, SP
                rr(i) = r(ilow) + DBLE(i-1)*deltar
                CALL CubicSplineInterpolation(r, U_r, U_r2, ilow, Nr, &
                                             &rr(i), UU(i))

             END DO

             !Perform HUNT again, finding ilow2 and ilow2+1 which are 
             !the indices that bound
             !where TargetedSValue is located
             flag = 0
             DO i=1, SP
                IF( UU(i) .GE. TargetedSValue ) THEN
                   ilow2 = i-1
                   flag = 1
                   EXIT
                END IF
             END DO

             IF( flag .EQ. 0) THEN           
                WRITE(*,*) 'Error in hunting for special points'
                !READ(*, '()')
                STOP
             END IF

             !Find the value of r where TargetedSValue is located 
             !by performing linear interpolation
             USp(ss) = rr(ilow2) + &
                  &( rr(ilow2+1) - rr(ilow2) )/( UU(ilow2+1) - UU(ilow2) ) *&
                  &(TargetedSValue - UU(ilow2) )

       END DO
       !$OMP END PARALLEL DO 
    
       RETURN
     END SUBROUTINE FindU
