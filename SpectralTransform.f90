!--------------------------------------------------------!
!   Subroutines for spectral transfomation               !
!--------------------------------------------------------!

    SUBROUTINE SpatialToSpectralTransform(&
            & Nr, Nth, Nphi,&
            & Mr, Mlm,&
            & rho, theta, phi,&
            & S, a)

        IMPLICIT    none
        INCLUDE     'shtns.f'

        !Calling variables
        INTEGER*4                   Nr, Nth, Nphi, Mr, Mlm
        REAL*8                      rho(Nr)
        REAL*8                      theta(Nth), phi(Nphi)
        COMPLEX*16                  S(Nr,Nth,Nphi)
        COMPLEX*16, INTENT(out) ::  a(Mr+1,Mlm)

        !Local variables
        INTEGER*4               i, j
        COMPLEX*16              f(Nr,Mlm)
    
        !Main subroutine
        
        !Use SHTns to find spectral coefficients at different r
        DO i=1,Nr
            CALL shtns_spat_to_sh(S(i,:,:), f(i,:))
        END DO
        !But now: f_lm(r) = Sum(over n) a_nlm T_n(r)
        !Transform to spectral coefficients a
        CALL ChebyshevSpatialToSpectral(Mlm,Nr,Mr,rho,f,a)

        RETURN
    END SUBROUTINE
!===========================================================!
    SUBROUTINE ChebyshevSpatialToSpectral(&
            & K, Nr, Mr,&
            & rho,&
            & f, a)

        !Using Chebyshev-Gauss quadrature assuming that rho is a grid
        !of Chebyshev roots

        USE         omp_lib
        IMPLICIT    none

        !Calling variables
        INTEGER*4                   K, Nr, Mr
        REAL*8                      rho(Nr)
        COMPLEX*16                  f(Nr,K)
        COMPLEX*16, INTENT(out)  :: a(Mr+1,K)
        
        !Local variables
        INTEGER*4               i,j,p,n
        REAL*8                  Spec2Spat(Nr, Mr+1)
        REAL*8                  wi

        !Main subroutine
        !Build wi (weight in the Chebyshev-Gauss quadrature*2/PI)
        wi = 2.0D0/DBLE(Mr+1)

        !$OMP PARALLEL 
        !$OMP DO PRIVATE(i)
        DO n=0,Mr
            DO i=1,Nr
                Spec2Spat(n+1,i) = wi * COS( n*ACOS(rho(i)) )
            END DO
        END DO
        !$OMP END DO

        !$OMP DO
        DO p=1,K
            a(:,p) = MATMUL(Spec2Spat, f(:,p))
        END DO
        !$OMP END DO
        !$OMP END PARALLEL

        RETURN
    END SUBROUTINE
!===========================================================!    
    SUBROUTINE SpectralToSpatialTransform(&
            & Nr, Nth, Nphi,&
            & Mr, Mlm,&
            & rho, theta, phi,&
            & a, S)

        IMPLICIT    none
        INCLUDE     'shtns.f'

        !Calling variables
        INTEGER*4                   Nr, Nth, Nphi, Mr, Mlm
        REAL*8                      rho(Nr)
        REAL*8                      theta(Nth), phi(Nphi)
        COMPLEX*16                  a(Mr+1,Mlm)
        COMPLEX*16, INTENT(OUT)  :: S(Nr,Nth,Nphi)

        !Local variables
        INTEGER*4               i, j
        COMPLEX*16              f(Nr,Mlm)
    
        !Main subroutine
        
        !f_lm(r) = Sum(over n) a_nlm T_n(r)
        CALL ChebyshevSpectralToSpatial(Mlm,Nr,Mr,rho,f,a)
    
        !Use SHTns to find get the spatial values
        DO i=1,Nr
            CALL shtns_sh_to_spat(f(i,:), S(i,:,:))
        END DO

        RETURN
    END SUBROUTINE
!===========================================================!    
    SUBROUTINE SpectralToSpatialOnARadialLine(&
            & Nr, Nth, Nphi,&
            & Mr, Mlm,&
            & rho, theta, phi,&
            & th0, phi0,&
            & a, S)

        IMPLICIT    none
        INCLUDE     'shtns.f'

        !Calling variables
        INTEGER*4                   Nr, Nth, Nphi, Mr, Mlm
        REAL*8                      rho(Nr)
        REAL*8                      theta(Nth), phi(Nphi)
        REAL*8                      th0, phi0
        COMPLEX*16                  a(Mr+1,Mlm)
        COMPLEX*16, INTENT(OUT)  :: S(Nr)

        !Local variables
        INTEGER*4               i, j
        COMPLEX*16              f(Nr,Mlm)
    
        !Main subroutine
        
        !f_lm(r) = Sum(over n) a_nlm T_n(r)
        CALL ChebyshevSpectralToSpatial(Mlm,Nr,Mr,rho,f,a)
    
        !Use SHTns to find get the spatial values
        DO i=1,Nr
            CALL shtns_sh_to_point(S(i), f(i,:), COS(th0), phi0)
        END DO

        RETURN
    END SUBROUTINE
!===========================================================!
    SUBROUTINE ChebyshevSpectralToSpatial(&
            & K, Nr, Mr,&
            & rho,&
            & f, a)

        !Using Chebyshev-Gauss quadrature assuming that rho is a grid
        !of Chebyshev roots

        USE         omp_lib
        IMPLICIT    none

        !Calling variables
        INTEGER*4                   K, Nr, Mr
        REAL*8                      rho(Nr)
        COMPLEX*16                  f(Nr,K)
        COMPLEX*16, INTENT(out)  :: a(Mr+1,K)
        
        !Local variables
        INTEGER*4               i,j,p,n
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
                
                Spec2Spat(i,n+1) = COS( n*ACOS(rho(i)) )

            END DO
        END DO
        !$OMP END DO

        !$OMP DO
        DO p=1,K
            f(:,p) = MATMUL(Spec2Spat, a(:,p))
        END DO
        !$OMP END DO
        !$OMP END PARALLEL

        RETURN
    END SUBROUTINE
!===========================================================!
    SUBROUTINE EvaluatedSdphi(&
            & Nr, Nth, Nphi,&
            & Mr, Mlm,&
            & rho, theta, phi,&
            & a, dSdphi)

        IMPLICIT     none 
        INCLUDE      'shtns.f'

        !Calling variables
        INTEGER*4                    Nr,Nth,Nphi,Mlm,Mr
        REAL*8                       rho(Nr), theta(Nth), phi(Nphi)
        COMPLEX*16                   a(Mr+1,Mlm)
        COMPLEX*16, INTENT(out)   :: dSdphi(Nr,Nth,Nphi)

        !Local variables
        INTEGER*4                    n, ml, l, lm
        COMPLEX*16                   ader(Mr+1,Mlm)

        !Main subroutine
        DO lm=1,Mlm
            CALL shtns_l_m(l,ml,lm)
            ader(:,lm) = CMPLX(0.0D0, 1.0D0)*DBLE(ml)*a(:,lm)
        END DO

        CALL SpectralToSpatialTransform(Nr,Nth,Nphi,Mr,Mlm,rho,theta,phi,&
                                       &ader,dSdphi)
        
        RETURN
    END SUBROUTINE
!===========================================================!
    SUBROUTINE EvaluatedSdr(&
            & Nr, Nth, Nphi,&
            & Mr, Mlm,&
            & rho, theta, phi,&
            & rmax, rmin,&
            & a, dSdr)

        IMPLICIT none

        !Calling variables
        INTEGER*4                    Nr,Nth,Nphi,Mlm,Mr
        REAL*8                       rmax,rmin
        REAL*8                       rho(Nr), theta(Nth), phi(Nphi)
        COMPLEX*16                   a(Mr+1,Mlm)
        COMPLEX*16, INTENT(out)   :: dSdr(Nr,Nth,Nphi)
        
        !Local variables
        REAL*8                       const, nder, nn
        COMPLEX*16                   ader(Mr+1,Mlm)

        !Main subroutine

        const = 2.0D0/(rmax-rmin)

        !Find derivatives by using recursion relations
        ader(Mr+1,:) = (0.0D0, 0.0D0)
        ader(Mr,:) = CMPLX(2.0D0*DBLE(Mr), 0.0D0)*a(Mr,:)
        DO nder = 1,(Mr-1)
            nn = Mr-nder+1
            ader(:,nn-1) = ader(:,nn+1) + &
                         & CMPLX(2.0D0*DBLE(nn-1),0.0D0)*a(nn,:)
        END DO
        ader = const*ader

        CALL SpectralToSpatialTransform(Nr,Nth,Nphi,Mr,Mlm,rho,theta,phi,&
                                       &ader,dSdr)

        RETURN
    END SUBROUTINE
!===========================================================!
    SUBROUTINE EvaluatedSdtheta(&
            & Nr, Nth, Nphi,&
            & Mr, Mlm,&
            & rho, theta, phi,&
            & a, dSdth)

        USE          omp_lib
        IMPLICIT     none 
        INCLUDE      'shtns.f'

        !Calling variables
        INTEGER*4                    Nr,Nth,Nphi,Mlm,Mr
        REAL*8                       rho(Nr), theta(Nth), phi(Nphi)
        COMPLEX*16                   a(Mr+1,Mlm)
        COMPLEX*16, INTENT(out)   :: dSdth(Nr,Nth,Nphi)

        !Local variables
        INTEGER*4                    i, j, k
        INTEGER*4                    n, ml, l, lm, mlp, lmp
        COMPLEX*16                   ader1(Mr+1,Mlm)
        COMPLEX*16                   ader2(Mr+1,Mlm)
        COMPLEX*16                   Term1(Nr,Nth,Nphi)
        COMPLEX*16                   Term2(Nr,Nth,Nphi)

        !Main subroutine

        ader2(:,:) = CMPLX(0.0D0, 0.0D0)
        
        DO lm=1,Mlm
            CALL shtns_l_m(l,ml,lm)
            ader1(:,lm) = DBLE(ml)*a(:,lm)
            
            IF( ml .LT. l ) THEN
                mlp = ml + 1
                CALL shtns_lmidx(lmp,l,mlp)
                ader2(:,lmp) = SQRT(DBLE((l-ml)*(l+ml+1))) * a(:,lm)
            END IF
        END DO

        CALL SpectralToSpatialTransform(Nr,Nth,Nphi,Mr,Mlm,rho,theta,phi,&
                                       &ader1,Term1)
        CALL SpectralToSpatialTransform(Nr,Nth,Nphi,Mr,Mlm,rho,theta,phi,&
                                       &ader2,Term2)

        !$OMP PARALLEL DO PRIVATE(j,k)
        DO i=1,Nr
          DO j=1,Nth
            DO k=1,Nphi
            
              dSdth(i,j,k) = Term1(i,j,k)/TAN(theta(j)) +&
                           & EXP(CMPLX(0.0D0, -1.0D0*phi(k)))*Term2(i,j,k)
            END DO
          END DO
        END DO
        !$OMP END PARALLEL DO

        RETURN
    END SUBROUTINE
