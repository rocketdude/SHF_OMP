!========================================================!
!    GetResults Subroutine                               !
!========================================================!

     SUBROUTINE GetResults(&
&Nr, Nth, Nphi, Mr, Lmax, Lgrid, SpM,&
&GLQWeights, GLQZeros,&
&r, rho, theta, phi,&
&a,&
&it, WriteSit,&
&thetaSp, phiSp,&
&TempAve)

!This subroutine finds U(theta, phi) using cubic spline interpolation
!and then performing linear interpolation

        USE             omp_lib
        IMPLICIT        none

!-----------------------------------------------------------!
!       Declare calling variables                           !
!-----------------------------------------------------------!

        INTEGER*4           :: Nr, Nth, Nphi, Mr, Lmax, Lgrid, SpM, it, WriteSit

        REAL*8              :: r(Nr), rho(Nr), theta(Nth), phi(Nphi)
        REAL*8              :: GLQWeights(Lgrid+1), GLQZeros(Lgrid+1)

        REAL*8              :: thetaSp(SpM)
        REAL*8              :: phiSp(SpM)
 
        COMPLEX*16          :: a(Mr+1,2,Lmax+1,Lmax+1)

        REAL*8, INTENT(out) :: TempAve

!--------------------------------------------------------!
!     Declare Locals                                     !
!--------------------------------------------------------!

        CHARACTER*32      CTemp

        INTEGER*4, PARAMETER ::  SP=1000 
        !Number of points to be calculated on the spline
        INTEGER*4, PARAMETER ::  PolyInterpOrder=6 
        !Number of points for poly. interp.
        INTEGER*4                i, j, k
        INTEGER*4                n, l, ml

        COMPLEX*16               S(Nr,Nth,Nphi)

!--------------------------------------------------------!
!      Main Subroutine                                   !
!--------------------------------------------------------!

        !Evaluate the eikonal data S
        CALL SpectralToSpatialTransform(Nr,Nth,Nphi,Mr,Lmax,Lgrid,&
                                      &GLQWeights,GLQZeros,&
                                      &rho,theta,phi,&
                                      &a,S)

        !Writing S into file (only at certain iterations)
        IF( MOD(it, WriteSit) .EQ. 0 ) CALL WriteS(Nr, Nth, Nphi, ABS(S), it)

        TempAve = 0.0D0
        DO i=1,Nr
            DO j=1,Nth
                DO k=1,Nphi
                   TempAve = TempAve + ABS( S(i,j,k) )
                END DO
            END DO
        END DO
        TempAve = TempAve / (Nr*Nth*Nphi)            
    
        RETURN
    END SUBROUTINE GetResults
