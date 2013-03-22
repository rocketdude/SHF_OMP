!========================================================!
!    Area Subroutine                                     !
!========================================================!

    SUBROUTINE CalculateArea(&
                &Nr,Nth,Nphi,&
                &Lmax,Lgrid,&
                &filterO,&
                &rho,theta,phi,&
                &rmax,rmin,&
                &GLQZeros,GLQWeights,&
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

    INTEGER*4               Nr,Nth,Nphi,Lmax,Lgrid,filterO
    REAL*8                  rho(Nr),theta(Nth),phi(Nphi),rmax,rmin
    REAL*8                  GLQZeros(Lgrid+1),GLQWeights(Lgrid+1)
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

    INTEGER*4               Nmax
    PARAMETER               (Nmax=6)

    COMPLEX*16              ulm(2,Lmax+1,Lmax+1)
    COMPLEX*16              dUdth(Nth,Nphi)
    COMPLEX*16              dUdphi(Nth,Nphi)
    REAL*8                  ththTerm,thphiTerm,phiphiTerm
    REAL*8                  r(Nr), g2D(Nr)
    REAL*8                  gu11,gu22,gu33
    REAL*8                  gu12,gu13,gu23
    REAL*8                  gd11,gd22,gd33
    REAL*8                  gd12,gd13,gd23
    REAL*8                  detGFunc(Nth,Nphi)
    REAL*8                  alm(2,Lmax+1,Lmax+1)
    REAL*8                  blm(2,Lmax+1,Lmax+1)
    REAL*8                  gN(Nr)
    REAL*8                  Urho
    REAL*8                  rdUdth,rdUdphi

    REAL*8                  G(Nmax)
    REAL*8                  rr(Nmax)

    INTEGER*4               i,j,k,il,iu
    INTEGER*4               error

!-----------------------------------------------------------!
!       Main Subroutine                                     !
!-----------------------------------------------------------!

    CALL GetRadialCoordinates(Nr,rmax,rmin,rho,r)

    CALL AngularToSpectralTransform(Nth,Nphi,Lmax,Lgrid,&
                &GLQWeights,GLQZeros,theta,phi,DCMPLX(U(:,:),0.0D0),ulm)
    CALL Evaluate2DdSdtheta(Nth,Nphi,Lmax,Lgrid,&
                            &GLQWeights,GLQZeros,filterO,theta,phi,ulm,dUdth)
    CALL Evaluate2DdSdphi(Nth,Nphi,Lmax,Lgrid,&
                            &GLQWeights,GLQZeros,filterO,theta,phi,ulm,dUdphi)

    !$OMP PARALLEL DO &
    !$OMP & PRIVATE(j,k,i,il,iu,g2D,gu11,gu22,gu33,gu12,gu13,gu23,&
    !$OMP & gd11,gd22,gd33,gd12,gd13,gd23,&
    !$OMP & G,rr,gN,Urho,rdUdth,rdUdphi,&
    !$OMP & ththTerm, thphiTerm, phiphiTerm)
    DO j=1,Nth
        DO k=1,Nphi
 
!!!            !Perform cubic spline interpolation
!!!            CALL ComputeSpline2ndDeriv(r,gRR(:,j,k),Nr,1.0D31,1.0D31,g2D)
!!!            CALL CubicSplineInterpolation(r,gRR(:,j,k),g2D,i,Nr,U(j,k),gu11)
!!!
!!!            CALL ComputeSpline2ndDeriv(r,gThTh(:,j,k),Nr,1.0D31,1.0D31,g2D)
!!!            CALL CubicSplineInterpolation(r,gThTh(:,j,k),g2D,i,Nr,U(j,k),gu22)
!!!
!!!            CALL ComputeSpline2ndDeriv(r,gPhiPhi(:,j,k),Nr,1.0D31,1.0D31,g2D)
!!!            CALL CubicSplineInterpolation(r,gPhiPhi(:,j,k),g2D,i,Nr,&
!!!                                        &U(j,k),gu33)
!!!
!!!            CALL ComputeSpline2ndDeriv(r,gRTh(:,j,k),Nr,1.0D31,1.0D31,g2D)
!!!            CALL CubicSplineInterpolation(r,gRTh(:,j,k),g2D,i,Nr,U(j,k),gu12)
!!!
!!!            CALL ComputeSpline2ndDeriv(r,gRPhi(:,j,k),Nr,1.0D31,1.0D31,g2D)
!!!            CALL CubicSplineInterpolation(r,gRPhi(:,j,k),g2D,i,Nr,U(j,k),gu13)
!!!
!!!            CALL ComputeSpline2ndDeriv(r,gThPhi(:,j,k),Nr,1.0D31,1.0D31,g2D)
!!!            CALL CubicSplineInterpolation(r,gThPhi(:,j,k),g2D,i,Nr,U(j,k),gu23)


!!!            !First we have to find where the U is located
!!!            CALL Locate(r,Nr,U(j,k),i)
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

            CALL GetRho(1,rmax,rmin,U(j,k),Urho)

            CALL SingleLineSpatialToSpectral(Nr,Nr-1,rho,gRR(:,j,k),gN)
            CALL SinglePointSpectralToSpatial(Nr-1,Urho,gN,gu11)

            CALL SingleLineSpatialToSpectral(Nr,Nr-1,rho,gThTh(:,j,k),gN)
            CALL SinglePointSpectralToSpatial(Nr-1,Urho,gN,gu22)

            CALL SingleLineSpatialToSpectral(Nr,Nr-1,rho,gPhiPhi(:,j,k),gN)
            CALL SinglePointSpectralToSpatial(Nr-1,Urho,gN,gu33)

            CALL SingleLineSpatialToSpectral(Nr,Nr-1,rho,gRTh(:,j,k),gN)
            CALL SinglePointSpectralToSpatial(Nr-1,Urho,gN,gu12)

            CALL SingleLineSpatialToSpectral(Nr,Nr-1,rho,gRPhi(:,j,k),gN)
            CALL SinglePointSpectralToSpatial(Nr-1,Urho,gN,gu13)

            CALL SingleLineSpatialToSpectral(Nr,Nr-1,rho,gThPhi(:,j,k),gN)
            CALL SinglePointSpectralToSpatial(Nr-1,Urho,gN,gu23)

            !Invert to get the lower index 3-metric
            CALL Invert3Metric(gu11,gu22,gu33,gu12,gu13,gu23,&
                              &gd11,gd22,gd33,gd12,gd13,gd23,&
                              &error)
            rdUdth = DBLE(dUdth(j,k))
            rdUdphi = DBLE(dUdphi(j,k))

            ththTerm = gd22 + gd11*rdUdth*rdUdth + 2.0D0*gd12*rdUdth
            thphiTerm = gd23 + gd11*rdUdth*rdUdphi +&
                      & gd12*rdUdphi + gd13*rdUdth
            phiphiTerm = gd33 + gd11*rdUdphi*rdUdphi +&
                      & 2.0D0*gd13*rdUdphi

            !detGFunc = sqrt( det(g) )/sin(theta) for theta,phi
            detGFunc(j,k) = &
                  &SQRT(ABS(ththTerm*phiphiTerm-thphiTerm*thphiTerm))&
                            &/SIN(theta(j))

        END DO
    END DO
    !$OMP END PARALLEL DO

    !Calculate area with Gauss-Legendre quadrature (FFT)
    CALL RealAngularToSpectralTransform(Nth,Nphi,Lmax,Lgrid,&
                    &GLQRealW,GLQRealZ,theta,phi,detGFunc,blm)
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

!========================================================!
!    Bisection Subroutine                                !
!========================================================!

    SUBROUTINE Bisect(yy,xx,n,xresult)

    IMPLICIT none

    INTEGER*4             n
    REAL*8                xx(n),yy(n)
    REAL*8,INTENT(out) :: xresult

    INTEGER*4,PARAMETER :: maxit = 200
    INTEGER*4,PARAMETER :: intOrder = 6 !must be even
    REAL*8,PARAMETER :: tol = 1.0D-16 !tolerance
    REAL*8,PARAMETER :: tol2 = 1.0D-14 !tolerance

    INTEGER*4 j,ilow
    REAL*8    dx,f,fmid,xmid,x1,x2
    REAL*8    yInt(intOrder),xInt(intOrder)

    CALL Locate(yy,n,0.0D0,ilow)
    IF((ilow .EQ. 0) .OR. (ilow .EQ. n)) STOP "***ERROR*** in bisect"

    fmid = yy(ilow+1)
    f = yy(ilow)

    x1 = xx(ilow)
    x2 = xx(ilow+1)

    IF(f*fmid .GE. 0.0D0) STOP "***ERROR*** root must be bracketed"
    
    IF(f .LT. 0.0D0) THEN
        xresult = x1
        dx = x2-x1
    ELSE
        xresult = x2
        dx = x1-x2
    END IF

    yInt(:) = yy( (ilow-(intOrder/2)-1):(ilow+(intOrder/2)) )
    xInt(:) = xx( (ilow-(intOrder/2)-1):(ilow+(intOrder/2)) )

    DO j=1,maxit
        dx = 0.5D0*dx
        xmid = xresult+dx

        !do 6th order interpolation
        CALL PolynomialInterpolation(xInt,yInt,intOrder,xmid,fmid)

        IF(fmid .LE. 0.0D0) xresult = xmid
        IF(ABS(dx).LT.tol .OR. fmid .EQ. 0.0D0) RETURN
    END DO

    IF(ABS(fmid) .GT. tol2) STOP "Bisection is too coarse"

    RETURN
    END SUBROUTINE Bisect
!========================================================!
!    Bisection Subroutine                                !
!========================================================!

    SUBROUTINE Bisect2(yy,xx,n,xresult)

    IMPLICIT none

    INTEGER*4             n
    REAL*8                xx(n),yy(n)
    REAL*8,INTENT(out) :: xresult

    INTEGER*4,PARAMETER :: maxit = 200
    REAL*8,PARAMETER :: tol = 1.0D-16 !tolerance
    REAL*8,PARAMETER :: tol2 = 1.0D-14 !tolerance

    INTEGER*4 j,ilow
    REAL*8    dx,f,fmid,xmid,x1,x2
    REAL*8    a(n)

    CALL Locate(yy,n,0.0D0,ilow)
    IF((ilow .EQ. 0) .OR. (ilow .EQ. n)) STOP "***ERROR*** in bisect"

    fmid = yy(ilow+1)
    f = yy(ilow)

    x1 = xx(ilow)
    x2 = xx(ilow+1)

    IF(f*fmid .GE. 0.0D0) STOP "***ERROR*** root must be bracketed"
    
    IF(f .LT. 0.0D0) THEN
        xresult = x1
        dx = x2-x1
    ELSE
        xresult = x2
        dx = x1-x2
    END IF

    CALL SingleLineSpatialToSpectral(n,n-1,xx,yy,a)

    DO j=1,maxit
        dx = 0.5D0*dx
        xmid = xresult+dx

        CALL SinglePointSpectralToSpatial(n-1,xmid,a,fmid)

        IF(fmid .LE. 0.0D0) xresult = xmid
        IF(ABS(dx).LT.tol .OR. fmid .EQ. 0.0D0) RETURN
    END DO

    IF(ABS(fmid) .GT. tol2) STOP "Bisection is too coarse"

    RETURN
    END SUBROUTINE Bisect2
