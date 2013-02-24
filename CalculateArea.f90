!========================================================!
!    Area Subroutine                                     !
!========================================================!

    SUBROUTINE CalculateArea(&
                &Nr,Nth,Nphi,&
                &Lmax,Lgrid,&
                &rho,theta,phi,&
                &rmax,rmin,&
                &GLQRealZ,GLQRealW,&
                &gRR,gThTh,gPhiPhi,&
                &gRTh,gRPhi,gThPhi,&
                &U,Area)

    !Calculates the area using Legendre-Gauss Quadrature.
    !The quadrature is done by using SHTOOLS

    USE         omp_lib
    USE         SHTOOLS
    IMPLICIT    none

!-----------------------------------------------------------!
!       Declare calling variables                           !
!-----------------------------------------------------------!

    INTEGER*4               Nr,Nth,Nphi,Lmax,Lgrid
    REAL*8                  rho(Nr),theta(Nth),phi(Nphi),rmax,rmin
    REAL*8                  GLQRealZ(Lgrid+1),GLQRealW(Lgrid+1)
    REAL*8                  gRR(Nr,Nth,Nphi)
    REAL*8                  gThTh(Nr,Nth,Nphi)
    REAL*8                  gPhiPhi(Nr,Nth,Nphi)
    REAL*8                  gRTh(Nr,Nth,Nphi)
    REAL*8                  gRPhi(Nr,Nth,Nphi)
    REAL*8                  gThPhi(Nr,Nth,Nphi)

    REAL*8                  U(Nth,Nphi)
    REAL*8,INTENT(out)  ::  Area

!-----------------------------------------------------------!
!       Declare local variables                             !
!-----------------------------------------------------------!

    REAL*8                  PI
    PARAMETER               (PI=3.141592653589793238462643383279502884197D0)
    REAL*8                  r(Nr), g2D(Nr)
    REAL*8                  gu11,gu22,gu33
    REAL*8                  gu12,gu13,gu23
    REAL*8                  gd11,gd22,gd33
    REAL*8                  gd12,gd13,gd23
    REAL*8                  detGFunc(Nth,Nphi)
    REAL*8                  blm(2,Lmax+1,Lmax+1)

    INTEGER*4               i,j,k,il,iu
    INTEGER*4               error

!-----------------------------------------------------------!
!       Main Subroutine                                     !
!-----------------------------------------------------------!

    CALL GetRadialCoordinates(Nr,rmax,rmin,rho,r)

    !$OMP PARALLEL DO &
    !$OMP & PRIVATE(j,k,i,il,iu,g2D,gu11,gu22,gu33,gu12,gu13,gu23,&
    !$OMP & gd11,gd22,gd33,gd12,gd13,gd23)
    DO j=1,Nth
        DO k=1,Nphi

            !First we have to find where the U is located
            CALL Locate(r,Nr,U(j,k),i)
 
            !Perform cubic spline interpolation
            CALL ComputeSpline2ndDeriv(r,gRR(:,j,k),Nr,1.0D31,1.0D31,g2D)
            CALL CubicSplineInterpolation(r,gRR(:,j,k),g2D,i,Nr,U(j,k),gu11)

            CALL ComputeSpline2ndDeriv(r,gThTh(:,j,k),Nr,1.0D31,1.0D31,g2D)
            CALL CubicSplineInterpolation(r,gThTh(:,j,k),g2D,i,Nr,U(j,k),gu22)

            CALL ComputeSpline2ndDeriv(r,gPhiPhi(:,j,k),Nr,1.0D31,1.0D31,g2D)
            CALL CubicSplineInterpolation(r,gPhiPhi(:,j,k),g2D,i,Nr,&
                                        &U(j,k),gu33)

            CALL ComputeSpline2ndDeriv(r,gRTh(:,j,k),Nr,1.0D31,1.0D31,g2D)
            CALL CubicSplineInterpolation(r,gRTh(:,j,k),g2D,i,Nr,U(j,k),gu12)

            CALL ComputeSpline2ndDeriv(r,gRPhi(:,j,k),Nr,1.0D31,1.0D31,g2D)
            CALL CubicSplineInterpolation(r,gRPhi(:,j,k),g2D,i,Nr,U(j,k),gu13)

            CALL ComputeSpline2ndDeriv(r,gThPhi(:,j,k),Nr,1.0D31,1.0D31,g2D)
            CALL CubicSplineInterpolation(r,gThPhi(:,j,k),g2D,i,Nr,U(j,k),gu23)

!!!            IF((i.EQ.0) .OR. (i.EQ.Nr)) THEN
!!!                PRINT *, 'Cannot find U in between r'
!!!                PRINT *, 'i =', i, ' and U(',j,',k',k,') =', U(j,k)
!!!                PRINT *, 'rmax = ', r(Nr), ',rmin = ', r(1)
!!!                STOP
!!!            ELSE IF(i .LE. 2) THEN
!!!                il = 1
!!!            ELSE IF(i .GE. (Nr-2)) THEN
!!!                il = Nr-3
!!!            ELSE
!!!                il = i-1
!!!            END IF
!!!            iu = il+Nmax-1
!!!            rr(:) = r(il:iu)
!!!            !Now perform polynomial interpolation
!!!            G(:) = gRR(il:iu,j,k)
!!!            CALL PolynomialInterpolation(rr,G,Nmax,U(j,k),gu11)
!!!            G(:) = gThTh(il:iu,j,k)
!!!            CALL PolynomialInterpolation(rr,G,Nmax,U(j,k),gu22)
!!!            G(:) = gPhiPhi(il:iu,j,k)
!!!            CALL PolynomialInterpolation(rr,G,Nmax,U(j,k),gu33)
!!!            G(:) = gRTh(il:iu,j,k)
!!!            CALL PolynomialInterpolation(rr,G,Nmax,U(j,k),gu12)
!!!            G(:) = gRPhi(il:iu,j,k)
!!!            CALL PolynomialInterpolation(rr,G,Nmax,U(j,k),gu13)
!!!            G(:) = gThPhi(il:iu,j,k)
!!!            CALL PolynomialInterpolation(rr,G,Nmax,U(j,k),gu23)

            !Invert to get the lower index 3-metric
            CALL Invert3Metric(gu11,gu22,gu33,gu12,gu13,gu23,&
                              &gd11,gd22,gd33,gd12,gd13,gd23,&
                              &error)

            !detGFunc = sqrt( det(g) )/sin(theta) for theta,phi
            detGFunc(j,k) = SQRT(ABS(gd22*gd33-gd23*gd23))/SIN(theta(j))

        END DO
    END DO
    !$OMP END PARALLEL DO

    !Calculate area with Gauss-Legendre quadrature (FFT)
    CALL SHExpandGLQ(blm,Lgrid,detGFunc,GLQRealW,&
                    &ZERO=GLQRealZ,NORM=4,LMAX_CALC=Lmax)
    Area = blm(1,1,1)*SQRT(4.0D0*PI)

    RETURN
    END SUBROUTINE CalculateArea

!========================================================!
!    Locate Subroutine                                   !
!========================================================!

    SUBROUTINE Locate(xx,n,x,j)

    !Perform bisection to find j where x is between xx(j) and xx(j+1).
    !Similar to that found in numerical recipes, if j=0 or j=n,
    !x is out of range.
    IMPLICIT none

    INTEGER*4           n,j
    REAL*8              xx(n),x

    INTEGER*4           jl,jm,ju

    jl=0
    ju=n+1

    !While-loop
10  IF((ju-jl) .GT. 1) THEN
        jm = (ju+jl)/2
        IF((xx(n) .GE. xx(1)) .EQV. (x .GE. xx(jm))) THEN
            jl = jm
        ELSE
            ju = jm
        END IF
    GOTO 10
    END IF

    IF(x .EQ. xx(1)) THEN
        j=1
    ELSE IF(x .EQ. xx(n)) THEN
        j=n-1
    ELSE
        j=jl
    END IF

    RETURN

    END SUBROUTINE Locate
