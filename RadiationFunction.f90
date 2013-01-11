!---------------------------------------------------------!
!   Radiation Boundary Condition                          !
!---------------------------------------------------------!

    SUBROUTINE  RadiationFunction(&
                &Nth, Nphi, Lmax, Lgrid, GLQWeights, GLQZeros,&
                &R, theta, phi,&
                &a, Fn)

        USE         omp_lib
        IMPLICIT    none

        !Declare passed variables
        INTEGER*4                       Nth, Nphi, Lmax, Lgrid
        REAL*8                          GLQWeights(Lgrid+1), GLQZeros(Lgrid+1)
        REAL*8                          R
        REAL*8                          theta(Nth), phi(Nphi)
        REAL*8                          a(2,Lmax+1,Lmax+1)
        REAL*8, INTENT(OUT)         ::  Fn(2,Lmax+1,Lmax+1)

        !Declare local variables
        INTEGER*4                       l, ml, j
        REAL*8                          SolPhi, eps, alp, kappa, sigma
        REAL*8                          epsStar, alpStar
        REAL*8                          PI
        REAL*8                          b(2,Lmax+1,Lmax+1)
        REAL*8                          bder(2,Lmax+1,Lmax+1)
        REAL*8                          c(2,Lmax+1,Lmax+1)

        REAL*8                          emit(Nth,Nphi)
        REAL*8                          absorb(Nth,Nphi)

        !Parameters related to the problem
        SolPhi      = 1366.0D0              !Solar constant in W/m^2
        eps         = 0.90D0                !emissivity
        alp         = 0.45D0                !absorptivity
        kappa       = 2.0D2                 !heat conductivity in W/m-K

        !Stefan-Boltzmann constant
        sigma       = 5.670373D-8
        PI = 3.141592653589793238462643383279502884197D0

        !Main subroutine
        !$OMP PARALLEL DO
        DO l = 0, Lmax
            b(:,l+1,:) = a(:,l+1,:)*( R**l )
            bder(:,l+1,:) = a(:,l+1,:)*DBLE(l)*( R**(l-1) )
        END DO
        !$OMP END PARALLEL DO

        epsStar = eps/kappa
        alpStar = alp/kappa

        !Calculate emission term
        CALL SpectralToAngularTransform(Nth,Nphi,Lmax,Lgrid,&
            &GLQWeights,GLQZeros,theta,phi,b,emit)

        emit = epsStar * sigma * emit*emit*emit*emit

        CALL AngularToSpectralTransform(Nth,Nphi,Lmax,Lgrid,&
            &GLQWeights,GLQZeros,theta,phi,emit,b)

        !Calculate absorption term
        !$OMP PARALLEL DO
        DO j=1,Nth
            IF( theta(j) .LE. PI/2.0D0 ) THEN
                absorb(j,:) = -1.0D0*alpStar*COS(theta(j))*SolPhi
            ELSE
                absorb(j,:) = 0.0D0
            END IF
        END DO
        !$OMP END PARALLEL DO

        CALL AngularToSpectralTransform(Nth,Nphi,Lmax,Lgrid,&
            &GLQWeights,GLQZeros,theta,phi,absorb,c)

        Fn = bder + b + c

        RETURN
    END SUBROUTINE
