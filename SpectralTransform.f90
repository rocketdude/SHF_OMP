!--------------------------------------------------------!
!   Subroutines for spectral transfomation               !
!--------------------------------------------------------!

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
        REAL*8                      S(Nth,Nphi)
        REAL*8, INTENT(out) ::      a(2,Lmax+1,Lmax+1)

        !Main subroutine
        
        CALL SHExpandGLQ(a, Lgrid, &
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
        REAL*8                      a(2,Lmax+1,Lmax+1)
        REAL*8, INTENT(OUT) ::      S(Nth,Nphi)
 
        !Main subroutine
        
        CALL MakeGridGLQ(S, a, Lgrid, &
                &ZERO=GLQZeros, NORM=1, LMAX_CALC=Lmax)

        RETURN
    END SUBROUTINE
