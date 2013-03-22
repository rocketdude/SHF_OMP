!========================================================!
!    FindU Subroutine                                    !
!========================================================!

     SUBROUTINE FindU(&
&Nr, Nth, Nphi, Mr, Lmax, Lgrid, SpM,&
&GLQWeights, GLQZeros,&
&gRR, gThTh, gPhiPhi,&
&gRTh, gRPhi, gThPhi,&
&rmax, rmin,&
&rho, theta, phi,&
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

       INTEGER*4            :: Nr, Nth, Nphi, Mr, Lmax, Lgrid, SpM, it, WriteSit

       REAL*8               :: rmax,rmin
       REAL*8               :: rho(Nr), theta(Nth), phi(Nphi)
       REAL*8               :: GLQWeights(Lgrid+1), GLQZeros(Lgrid+1)

       REAL*8               :: thetaSp(SpM)
       REAL*8               :: phiSp(SpM)

       REAL*8               :: gRR(Nr,Nth,Nphi)
       REAL*8               :: gThTh(Nr,Nth,Nphi)
       REAL*8               :: gPhiPhi(Nr,Nth,Nphi)
       REAL*8               :: gRTh(Nr,Nth,Nphi)
       REAL*8               :: gRPhi(Nr,Nth,Nphi)
       REAL*8               :: gThPhi(Nr,Nth,Nphi)
    
       COMPLEX*16           :: a(Mr+1,2,Lmax+1,Lmax+1)

       REAL*8, INTENT(out)  :: U(Nth, Nphi)
       REAL*8, INTENT(out)  :: Uave
       REAL*8, INTENT(out)  :: USp(SpM)
       REAL*8, INTENT(out)  :: gRRUsqrd(Nth, Nphi)
       REAL*8, INTENT(out)  :: gRRUsqrdAve

!--------------------------------------------------------!
!     Declare Locals                                     !
!--------------------------------------------------------!

       CHARACTER*32      CTemp

       INTEGER*4, PARAMETER ::  Ncut=4
       REAL*8,PARAMETER ::    TargetedSValue=100.0D0
       INTEGER*4                i, j, k, ii, dii, jj, ss

       COMPLEX*16             S(Nr,Nth,Nphi)
    
       REAL*8                 r(Nr)
       REAL*8                 SLine(SpM,Nr)
       REAL*8                 U_r(Nr) 
       REAL*8                 U_r2(Nr-Ncut-Ncut)
       REAL*8                 r2(Nr-Ncut-Ncut)
 

!--------------------------------------------------------!
!      Main Subroutine                                   !
!--------------------------------------------------------!

       !not used
       gRRUSqrd = 0.0D0
       gRRUsqrdAve = 0.0D0

       !Evaluate the eikonal data S
       CALL SpectralToSpatialTransform(Nr,Nth,Nphi,Mr,Lmax,Lgrid,&
                                      &GLQWeights,GLQZeros,&
                                      &rho,theta,phi,&
                                      &a,S)

       !Writing S into file (only at certain iterations)
       IF( MOD(it, WriteSit) .EQ. 0 ) CALL WriteS(Nr, Nth, Nphi, DBLE(S),it)

       CALL GetRadialCoordinates(Nr,rmax,rmin,rho,r)

       !$OMP PARALLEL DO &
       !$OMP &PRIVATE(j,k,U_r,U_r2,r2)
       DO j = 1, Nth
          DO k = 1, Nphi

            U_r = DBLE(S(:,j,k)) - TargetedSValue
            U_r2(:) = U_r((1+Ncut):(Nr-Ncut))
            r2(:) = r((1+Ncut):(Nr-Ncut))

            !Now perform bisection
            !CALL Bisect(U_r2,r2,(Nr-2*Ncut),U(j,k))
            CALL Bisect2(U_r,rho,Nr,U(j,k))
            CALL GetRadialCoordinates(1,rmax,rmin,U(j,k),U(j,k))

          END DO
       END DO
       !$OMP END PARALLEL DO

       !Calculate the average U: Uave 
       !(only applies if there's obvious spherical symmetry)
       Uave = SUM( U )/ (Nth*Nphi)

       CALL GetRealSpatialValueOnRadialLine(Nr,Nth,Nphi,Mr,SpM,Lmax,Lgrid,&
                                           &rho,theta,phi,thetaSp,phiSp,&
                                           &S, SLine)

       !Calculate the values of U at the requested directions
       !$OMP PARALLEL DO &
       !$OMP &PRIVATE(ss,U_r,U_r2,r2)
       DO ss = 1, SpM
            U_r = SLine(ss,:) - TargetedSValue
            U_r2(:) = U_r((1+Ncut):(Nr-Ncut))
            r2(:) = r((1+Ncut):(Nr-Ncut))

            !Now perform bisection
            !CALL Bisect(U_r2,r2,(Nr-2*Ncut),USp(ss)) 
            CALL Bisect2(U_r,rho,Nr,USp(ss))
            CALL GetRadialCoordinates(1,rmax,rmin,USp(ss),USp(ss))
       END DO
       !$OMP END PARALLEL DO 
      
       RETURN
     END SUBROUTINE FindU
