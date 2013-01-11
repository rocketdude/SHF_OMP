!--------------------------------------------------------!
!    GetResults subroutine                               !
!--------------------------------------------------------!

      SUBROUTINE GetResults(& 
& Nth, Nphi,&
& Lmax, Lgrid,&
& GLQWeights, GLQZeros,&
& R, theta, phi,&
& a, S)

        !This subroutine calculates the initial eikonal data S

        USE            omp_lib
        IMPLICIT       none

!--------------------------------------------------------!
!     Declare calling variables                          !
!--------------------------------------------------------!

        INTEGER*4               :: Nth, Nphi
        INTEGER*4               :: Lmax, Lgrid
        
        REAL*8                  :: GLQWeights(Lgrid+1), GLQZeros(Lgrid+1)
        REAL*8                  :: R
        REAL*8                  :: theta(Nth)
        REAL*8                  :: phi(Nphi)

        REAL*8                  :: a(2, Lmax+1, Lmax+1)
        REAL*8, INTENT(OUT)     :: S(Nth,Nphi)

!--------------------------------------------------------!
!     Declare Locals                                     !
!--------------------------------------------------------!

        INTEGER*4               l

!--------------------------------------------------------!
!      Main Subroutine                                   !
!--------------------------------------------------------!

        !$OMP PARALLEL DO
        DO l=0,Lmax
            a(:,l+1,:) = a(:,l+1,:) * (R**l)
        END DO
        !$OMP END PARALLEL DO
        
        CALL SpectralToAngularTransform(Nth,Nphi,Lmax,Lgrid,&
                                    &GLQWeights,GLQZeros,theta,phi,a,S)

        PRINT *, 'Maximum temperature =', MAXVAL(S)
        PRINT *, 'Minimum temperature =', MINVAL(S)

        RETURN
      END SUBROUTINE GetResults
