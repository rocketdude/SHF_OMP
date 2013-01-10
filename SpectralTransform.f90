!--------------------------------------------------------!
!   Subroutines for spectral transfomation               !
!--------------------------------------------------------!

    SUBROUTINE SpatialToSpectralTransform(&
            & Nr, Nth, Nphi,&
            & Mr, Lmax, Lgrid,&
            & GLQWeights, GLQZeros,&
            & rho, theta, phi,&
            & S, a)

        USE SHTOOLS
        IMPLICIT    none

        !Calling variables
        INTEGER*4                   Nr, Nth, Nphi, Mr, Lmax, Lgrid
        REAL*8                      rho(Nr)
        REAL*8                      theta(Nth), phi(Nphi)
        REAL*8                      GLQWeights(Lgrid+1), GLQZeros(Lgrid+1)
        COMPLEX*16                  S(Nr,Nth,Nphi)
        COMPLEX*16, INTENT(out) ::  a(Mr+1,2,Lmax+1,Lmax+1)

        !Local variables
        INTEGER*4               i, j
        COMPLEX*16              f(Nr,2,Lmax+1,Lmax+1)
    
        !Main subroutine
        
        !Use SHTOOLS to find spectral coefficients at different r
        DO i=1,Nr
            CALL SHExpandGLQC(f(i,:,:,:), Lgrid, &
                    &S(i,:,:), GLQWeights, ZERO=GLQZeros,NORM=1,LMAX_CALC=Lmax)
        END DO
        !But now: f_lm(r) = Sum(over n) a_nlm T_n(r)
        !Transform to spectral coefficients a
        CALL ChebyshevSpatialToSpectral(Lmax+1,Nr,Mr,rho,f,a)

        RETURN
    END SUBROUTINE
!===========================================================!
    SUBROUTINE ChebyshevSpatialToSpectral(&
            & K, Nr, Mr,&
            & rho,&
            & f, a)

        !Using Chebyshev-Gauss-Radau quadrature
        !Note that rho and r are inverted

        USE         omp_lib
        IMPLICIT    none

        !Calling variables
        INTEGER*4                   K, Nr, Mr
        REAL*8                      rho(Nr)
        COMPLEX*16                  f(Nr,2,K,K)
        COMPLEX*16, INTENT(out)  :: a(Mr+1,2,K,K)
        
        !Local variables
        INTEGER*4               i,j,p,q,n
        REAL*8                  Spec2Spat(Nr, Mr+1)
        REAL*8                  wi

        !Main subroutine

        !$OMP PARALLEL 
        !$OMP DO PRIVATE(i, wi)
        DO n=0,Mr
            DO i=1,Nr
                
                !wi is the weight in Chebyshev-Gauss-Radau quadrature
                IF( i .EQ. 1 ) THEN
                    wi = 2.0D0 / ( 2.0D0 * DBLE(Nr-1) + 1.0D0 )
                ELSE 
                    wi = 4.0D0 / ( 2.0D0 * DBLE(Nr-1) + 1.0D0 )
                END IF

                Spec2Spat(n+1,i) = wi*COS( DBLE(n)*ACOS(rho(i)) )
            END DO
        END DO
        !$OMP END DO

        !$OMP DO
        DO p=1,K
            DO q=1,K
                a(:,1,p,q) = MATMUL( Spec2Spat, f(:,1,p,q) )
                a(:,2,p,q) = MATMUL( Spec2Spat, f(:,2,p,q) )
            END DO
        END DO
        !$OMP END DO
        !$OMP END PARALLEL

        RETURN
    END SUBROUTINE
!===========================================================!    
    SUBROUTINE SpectralToSpatialTransform(&
            & Nr, Nth, Nphi,&
            & Mr, Lmax, Lgrid,&
            & GLQWeights, GLQZeros,&            
            & rho, theta, phi,&
            & a, S)

        USE SHTOOLS
        IMPLICIT    none

        !Calling variables
        INTEGER*4                   Nr, Nth, Nphi, Mr, Lmax, Lgrid
        REAL*8                      GLQWeights(Lgrid+1), GLQZeros(Lgrid+1)
        REAL*8                      rho(Nr)
        REAL*8                      theta(Nth), phi(Nphi)
        COMPLEX*16                  a(Mr+1,2,Lmax+1,Lmax+1)
        COMPLEX*16, INTENT(OUT)  :: S(Nr,Nth,Nphi)

        !Local variables
        INTEGER*4               i, j
        COMPLEX*16              f(Nr,2,Lmax+1,Lmax+1)
    
        !Main subroutine
        
        !f_lm(r) = Sum(over n) a_nlm T_n(r)
        CALL ChebyshevSpectralToSpatial(Lmax+1,Nr,Mr,rho,a,f)
    
        !Use SHTOOLS to find get the spatial values
        DO i=1,Nr
            CALL MakeGridGLQC(S(i,:,:), f(i,:,:,:), Lgrid, &
                    &ZERO=GLQZeros, NORM=1, LMAX_CALC=Lmax)
        END DO

        RETURN
    END SUBROUTINE
!===========================================================!    
    SUBROUTINE GetRealSpatialValueOnRadialLine(&
            & Nr, Nth, Nphi,&
            & Mr, SpM, Lmax, Lgrid,&
            & rho, theta, phi,&
            & th0, phi0,&
            & gridS, S)

        USE SHTOOLS
        IMPLICIT    none

        !Calling variables
        INTEGER*4                   Nr, Nth, Nphi, Mr, Lmax, Lgrid, SpM
        REAL*8                      rho(Nr)
        REAL*8                      theta(Nth), phi(Nphi)
        REAL*8                      th0(SpM), phi0(SpM)
        COMPLEX*16                  gridS(Nr,Nth,Nphi)
        REAL*8, INTENT(OUT)  ::     S(SpM,Nr)

        !Local variables
        INTEGER*4               i, j
        REAL*8                  f(Nr,2,Lmax+1,Lmax+1)
        REAL*8                  SReal(Nr,Nth,Nphi)
        REAL*8                  Zeros(Lgrid+1), W(Lgrid+1)
        REAL*8                  th0Degrees(SpM), phi0Degrees(SpM)
        REAL*8, PARAMETER ::    PI = 3.141592653589793238462643383279502884197D0
    
        !Main subroutine
         
        !SHTOOLS only accept degrees
        th0Degrees(:) = 90.0D0 - (th0(:)*180.0D0/PI)
        phi0Degrees(:) = phi0(:) * 180.0D0 / PI
        SReal = ABS(gridS)

        CALL PreCompute(Lgrid, Zeros, W)
        
        !Get the real spherical harmonic coefficients
        DO i=1,Nr
            CALL SHExpandGLQ(f(i,:,:,:), Lgrid, SReal(i,:,:), W, &
                    &ZERO=Zeros, NORM=1, LMAX_CALC=Lmax)
        END DO

        DO j=1,SpM
            DO i=1,Nr
                !Use these coefficients f to find the value of S 
                S(j,i) = MakeGridPoint(f(i,:,:,:), Lmax, &
                        & th0Degrees(j), phi0Degrees(j), NORM=1 )
            END DO
        END DO

        RETURN
    END SUBROUTINE
!===========================================================!
    SUBROUTINE ChebyshevSpectralToSpatial(&
            & K, Nr, Mr,&
            & rho,&
            & a, f)

        !Using Chebyshev-Gauss-Radau Quadrature
        !Note that rho and r are inverted

        USE         omp_lib
        IMPLICIT    none

        !Calling variables
        INTEGER*4                   K, Nr, Mr
        REAL*8                      rho(Nr)
        COMPLEX*16                  a(Mr+1,2,K,K)
        COMPLEX*16, INTENT(out)  :: f(Nr,2,K,K)
        
        !Local variables
        INTEGER*4               i,j,p,q,n
        REAL*8                  const
        REAL*8                  Spec2Spat(Mr+1, Nr)

        !Main subroutine

        !$OMP PARALLEL 
        !$OMP DO PRIVATE(i)
        DO i=1,Nr
            DO n=0,Mr
                IF(n.EQ.0) THEN
                    const = 0.5D0
                ELSE
                    const = 1.0D0
                END IF
                
                Spec2Spat(i,n+1) = const*COS( DBLE(n)*ACOS(rho(i)) )

            END DO
        END DO
        !$OMP END DO

        !$OMP DO
        DO p=1,K
            DO q=1,K
                f(:,1,p,q) = MATMUL(Spec2Spat, a(:,1,p,q))
                f(:,2,p,q) = MATMUL(Spec2Spat, a(:,2,p,q))
            END DO
        END DO
        !$OMP END DO
        !$OMP END PARALLEL

        RETURN
    END SUBROUTINE
!===========================================================!
    SUBROUTINE EvaluatedSdphi(&
            & Nr, Nth, Nphi,&
            & Mr, Lmax, Lgrid,&
            & GLQWeights, GLQZeros,&            
            & rho, theta, phi,&
            & a, dSdphi)

        USE SHTOOLS
        IMPLICIT     none 

        !Calling variables
        INTEGER*4                    Nr,Nth,Nphi,Mr,Lmax,Lgrid
        REAL*8                       GLQWeights(Lgrid+1),GLQZeros(Lgrid+1)
        REAL*8                       rho(Nr), theta(Nth), phi(Nphi)
        COMPLEX*16                   a(Mr+1,2,Lmax+1,Lmax+1)
        COMPLEX*16, INTENT(out)   :: dSdphi(Nr,Nth,Nphi)

        !Local variables
        INTEGER*4                    n,l,ml
        COMPLEX*16                   ader(Mr+1,2,Lmax+1,Lmax+1)

        ader(:,:,:,:) = DCMPLX(0.0D0, 0.0D0)

        !Main subroutine
        DO l=0,Lmax
            DO ml=0,l
            
            IF( ml .EQ. 0 ) THEN
                ader(:,1,l+1,ml+1) = DCMPLX(0.0D0, 1.0D0)*DBLE(ml)&
                                    &*a(:,1,l+1,ml+1)
            ELSE
                ader(:,1,l+1,ml+1) = DCMPLX(0.0D0, 1.0D0)*DBLE(ml)&
                                    &*a(:,1,l+1,ml+1)
                ader(:,2,l+1,ml+1) = DCMPLX(0.0D0,-1.0D0)*DBLE(ml)&
                                    &*a(:,2,l+1,ml+1)
            END IF
            
            END DO
        END DO

        CALL SpectralToSpatialTransform(Nr,Nth,Nphi,Mr,Lmax,Lgrid,&
                                       &GLQWeights,GLQZeros,&
                                       &rho,theta,phi,&
                                       &ader,dSdphi)
        
        RETURN
    END SUBROUTINE
!===========================================================!
    SUBROUTINE EvaluatedSdr(&
            & Nr, Nth, Nphi,&
            & Mr, Lmax, Lgrid,&
            & GLQWeights, GLQZeros,&            
            & rho, theta, phi,&
            & rmax, rmin,&
            & a, dSdr)

        IMPLICIT none

        !Calling variables
        INTEGER*4                    Nr,Nth,Nphi,Mr,Lmax,Lgrid
        REAL*8                       GLQWeights(Lgrid+1),GLQZeros(Lgrid+1)
        REAL*8                       rmax,rmin
        REAL*8                       rho(Nr), theta(Nth), phi(Nphi)
        COMPLEX*16                   a(Mr+1,2,Lmax+1,Lmax+1)
        COMPLEX*16, INTENT(out)   :: dSdr(Nr,Nth,Nphi)
        
        !Local variables
        REAL*8                       const
        INTEGER*4                    nder, nn
        COMPLEX*16                   ader(Mr+1,2,Lmax+1,Lmax+1)

        !Main subroutine
    
        !Note that r and rho are inverted because of the Chebyshev-Gauss-Radau
        const = -2.0D0/(rmax-rmin)

        !Find derivatives by using recursion relations
        !This is similar to chder from Numerical Recipes
        ader(Mr+1,:,:,:) = DCMPLX(0.0D0, 0.0D0)
        ader(Mr,:,:,:) = DCMPLX(2.0D0*DBLE(Mr), 0.0D0)*a(Mr,:,:,:)
        DO nder = 1,(Mr-1)
            nn = Mr-nder+1
            ader(nn-1,:,:,:) = ader(nn+1,:,:,:) + &
                         & DCMPLX(2.0D0*DBLE(nn-1),0.0D0)*a(nn,:,:,:)
        END DO
        ader = const*ader

        CALL SpectralToSpatialTransform(Nr,Nth,Nphi,Mr,Lmax,Lgrid,&
                                       &GLQWeights,GLQZeros,&
                                       &rho,theta,phi,&
                                       &ader,dSdr)

        RETURN
    END SUBROUTINE
!===========================================================!
    SUBROUTINE EvaluatedSdtheta(&
            & Nr, Nth, Nphi,&
            & Mr, Lmax, Lgrid,&
            & GLQWeights, GLQZeros,&            
            & rho, theta, phi,&
            & a, dSdth)

        USE          omp_lib
        USE          SHTOOLS
        IMPLICIT     none 

        !Calling variables
        INTEGER*4                    Nr,Nth,Nphi,Mr,Lmax,Lgrid
        REAL*8                       GLQWeights(Lgrid+1), GLQZeros(Lgrid+1)
        REAL*8                       rho(Nr), theta(Nth), phi(Nphi)
        COMPLEX*16                   a(Mr+1,2,Lmax+1,Lmax+1)
        COMPLEX*16, INTENT(out)   :: dSdth(Nr,Nth,Nphi)

        !Local variables
        INTEGER*4                    i, j, k
        INTEGER*4                    n, ml, l, mlp1, mlp2
        REAL*8                       const1, const2
        COMPLEX*16                   ader1(Mr+1,2,Lmax+1,Lmax+1)
        COMPLEX*16                   ader2(Mr+1,2,Lmax+1,Lmax+1)
        COMPLEX*16                   Term1(Nr,Nth,Nphi)
        COMPLEX*16                   Term2(Nr,Nth,Nphi)

        !Main subroutine

        ader1(:,:,:,:) = DCMPLX(0.0D0, 0.0D0)
        ader2(:,:,:,:) = DCMPLX(0.0D0, 0.0D0)
        
        DO l=0,Lmax
            DO ml=0,l

            IF( ml .EQ. 0 ) THEN
                ader1(:,1,l+1,ml+1) = DBLE(ml)*a(:,1,l+1,ml+1)
            ELSE
                ader1(:,1,l+1,ml+1) = DBLE(ml)*a(:,1,l+1,ml+1)
                ader1(:,2,l+1,ml+1) = -1.0D0*DBLE(ml)*a(:,2,l+1,ml+2)
            END IF
            
            const1 = SQRT(DBLE((l-ml)*(l+ml+1)))
            const2 = SQRT(DBLE((l+ml)*(l-ml+1)))
            mlp1 = ml + 1
            IF( ml .LT. l ) THEN
                ader2(:,1,l+1,mlp1+1) = const1*a(:,1,l+1,ml+1)
            END IF

            mlp2 = ABS(-1*ml + 1)
            IF( mlp2 .EQ. 0 ) THEN
                ader2(:,1,l+1,mlp2+1) = const2*a(:,2,l+1,ml+1)
            ELSE
                ader2(:,2,l+1,mlp2+1) = const2*a(:,2,l+1,ml+1)
            END IF

            END DO
        END DO

        CALL SpectralToSpatialTransform(Nr,Nth,Nphi,Mr,Lmax,Lgrid,&
                                       &GLQWeights,GLQZeros,&
                                       &rho,theta,phi,&
                                       &ader1,Term1)
        CALL SpectralToSpatialTransform(Nr,Nth,Nphi,Mr,Lmax,Lgrid,&
                                       &GLQWeights,GLQZeros,&
                                       &rho,theta,phi,&
                                       &ader2,Term2)

        !$OMP PARALLEL DO PRIVATE(j,k)
        DO i=1,Nr
          DO j=1,Nth
            DO k=1,Nphi
            
              dSdth(i,j,k) = Term1(i,j,k)/TAN(theta(j)) +&
                           & EXP(DCMPLX(0.0D0, -1.0D0*phi(k)))*Term2(i,j,k)
            END DO
          END DO
        END DO
        !$OMP END PARALLEL DO

        RETURN
    END SUBROUTINE
!===========================================================!
    SUBROUTINE Evaluated2Sdr2(&
            & Nr, Nth, Nphi,&
            & Mr, Lmax, Lgrid,&
            & GLQWeights, GLQZeros,&            
            & rho, theta, phi,&
            & rmax, rmin,&
            & dSdr, d2Sdr2)

        IMPLICIT none

        !Calling variables
        INTEGER*4                    Nr,Nth,Nphi,Mr,Lmax,Lgrid
        REAL*8                       GLQWeights(Lgrid+1),GLQZeros(Lgrid+1)
        REAL*8                       rmax,rmin
        REAL*8                       rho(Nr), theta(Nth), phi(Nphi)
        COMPLEX*16                   dSdr(Nr,Nth,Nphi)
        COMPLEX*16, INTENT(out)   :: d2Sdr2(Nr,Nth,Nphi)

        !Local variables
        COMPLEX*16                   dadr(Mr+1,2,Lmax+1,Lmax+1)
 
        CALL SpatialToSpectralTransform(Nr,Nth,Nphi,Mr,Lmax,Lgrid,&
                                       &GLQWeights,GLQZeros,&
                                       &rho,theta,phi,&
                                       &dSdr,dadr)
        CALL EvaluatedSdr(Nr,Nth,Nphi,Mr,Lmax,Lgrid,GLQWeights,GLQZeros,&
                &rho,theta,phi,rmax,rmin,dadr,d2Sdr2)

        RETURN
    END SUBROUTINE
!===========================================================!
    SUBROUTINE EvaluateAngularLaplacian(&
            & Nr, Nth, Nphi,&
            & Mr, Lmax, Lgrid,&
            & GLQWeights, GLQZeros,&            
            & rho, theta, phi,&
            & a, delsqrdS)

        USE SHTOOLS
        IMPLICIT     none 

        !Calling variables
        INTEGER*4                    Nr,Nth,Nphi,Mr,Lmax,Lgrid
        REAL*8                       GLQWeights(Lgrid+1),GLQZeros(Lgrid+1)
        REAL*8                       rho(Nr), theta(Nth), phi(Nphi)
        COMPLEX*16                   a(Mr+1,2,Lmax+1,Lmax+1)
        COMPLEX*16, INTENT(out)   :: delsqrdS(Nr,Nth,Nphi)

        !Local variables
        INTEGER*4                    n,l,ml
        COMPLEX*16                   ader(Mr+1,2,Lmax+1,Lmax+1)

        ader(:,:,:,:) = DCMPLX(0.0D0, 0.0D0)

        !Main subroutine
        DO l=0,Lmax
            DO ml=0,l
            
            IF( ml .EQ. 0 ) THEN
                ader(:,1,l+1,ml+1) = -1.0D0*DBLE(l*(l+1))*&
                                    &a(:,1,l+1,ml+1)
            ELSE
                ader(:,1,l+1,ml+1) = -1.0D0*DBLE(l*(l+1))*&
                                    &a(:,1,l+1,ml+1)
                ader(:,2,l+1,ml+1) = -1.0D0*DBLE(l*(l+1))*&
                                    &a(:,2,l+1,ml+1)
            END IF
            
            END DO
        END DO

        CALL SpectralToSpatialTransform(Nr,Nth,Nphi,Mr,Lmax,Lgrid,&
                                       &GLQWeights,GLQZeros,&
                                       &rho,theta,phi,&
                                       &ader,delsqrdS)
        
        RETURN
    END SUBROUTINE
!===========================================================!
    SUBROUTINE GetThetaAndPhi(Nth, Nphi, Lgrid, theta, phi)
    
        USE          SHTOOLS
        USE          omp_lib
        IMPLICIT     none 

        !Calling variables
        INTEGER*4                    Nth, Nphi, Lgrid
        REAL*8                       theta(Nth), phi(Nphi)
        REAL*8, PARAMETER ::    PI = 3.141592653589793238462643383279502884197D0

        !Local variables
        INTEGER*4                    NthOut, NphiOut, j,k
        REAL*8                       latitude(Nth), longitude(Nphi)
        
        !Main subroutine

        CALL GLQGridCoord(latitude,longitude,Lgrid,NthOut,NphiOut)
    
        IF( NthOut .NE. Nth ) STOP "***ERROR*** Nth is not correct"
        IF( NphiOut .NE. Nphi ) STOP "***ERROR*** Nphi is not correct"

        theta = ABS( latitude - 90D0 )*PI/180.0D0
        phi = longitude*PI/180.0D0

    END SUBROUTINE
!===========================================================!
    SUBROUTINE AngularToSpectralTransform(&
            & Nth, Nphi,&
            & Lmax, Lgrid,&
            & GLQWeights, GLQZeros,&
            & theta, phi,&
            & S, a)

        USE SHTOOLS
        IMPLICIT    none

        !Calling variables
        INTEGER*4                   Nth, Nphi, Lmax, Lgrid
        REAL*8                      theta(Nth), phi(Nphi)
        REAL*8                      GLQWeights(Lgrid+1), GLQZeros(Lgrid+1)
        COMPLEX*16                  S(Nth,Nphi)
        COMPLEX*16, INTENT(out) ::  a(2,Lmax+1,Lmax+1)

        !Main subroutine
        
        CALL SHExpandGLQC(a, Lgrid, &
                &S, GLQWeights, ZERO=GLQZeros,NORM=1,LMAX_CALC=Lmax)

        RETURN
    END SUBROUTINE
!===========================================================!
    SUBROUTINE SpectralToAngularTransform(&
            & Nth, Nphi,&
            & Lmax, Lgrid,&
            & GLQWeights, GLQZeros,&            
            & theta, phi,&
            & a, S)

        USE SHTOOLS
        IMPLICIT    none

        !Calling variables
        INTEGER*4                   Nth, Nphi, Lmax, Lgrid
        REAL*8                      GLQWeights(Lgrid+1), GLQZeros(Lgrid+1)
        REAL*8                      theta(Nth), phi(Nphi)
        COMPLEX*16                  a(2,Lmax+1,Lmax+1)
        COMPLEX*16, INTENT(OUT)  :: S(Nth,Nphi)
 
        !Main subroutine
         
        CALL MakeGridGLQC(S, a, Lgrid, &
                &ZERO=GLQZeros, NORM=1, LMAX_CALC=Lmax)

        RETURN
    END SUBROUTINE
