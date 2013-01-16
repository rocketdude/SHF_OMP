!--------------------------------------------------------!
!    GetInitialData subroutine                           !
!--------------------------------------------------------!

      SUBROUTINE GetInitialData(& 
& Nth, Nphi,&
& Lmax, Lgrid,&
& GLQWeights, GLQZeros,&
& T0, tolA, eps,&
& R, theta, phi,&
& a)

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
        REAL*8                  :: T0, tolA, eps

        REAL*8, INTENT(out)     :: a(2, Lmax+1, Lmax+1)

!--------------------------------------------------------!
!     Declare Locals                                     !
!--------------------------------------------------------!

        INTEGER*4        l, ml
        REAL*8           S(Nth,Nphi)


!--------------------------------------------------------!
!      Main Subroutine                                   !
!--------------------------------------------------------!

        S = T0

        CALL AngularToSpectralTransform(Nth,Nphi,Lmax,Lgrid,&
                                    &GLQWeights,GLQZeros,theta,phi,S,a)

        !$OMP PARALLEL DO
        DO l=0,Lmax
            a(:,l+1,:) = a(:,l+1,:) / (R**l) 
        END DO
        !$OMP END PARALLEL DO 

        !$OMP PARALLEL DO
        DO l=0,Lmax
            DO ml=0,Lmax
                IF( a(1,l+1,ml+1) .EQ. 0.0D0 ) a(1,l+1,ml+1) = tolA
                IF( a(2,l+1,ml+1) .EQ. 0.0D0 ) a(2,l+1,ml+1) = tolA
            END DO 
        END DO
        !$OMP END PARALLEL DO
        
        RETURN
      END SUBROUTINE GetInitialData
