!--------------------------------------------------------!
!    CalculateInitialData subroutine                     !
!--------------------------------------------------------!

      SUBROUTINE ReinitializeData(& 
& Nr, Nth, Nphi, &
& Mr, Lmax, Lgrid,&
& GLQWeights, GLQZeros,&
& rmax, rmin,&
& rho, theta, phi,&
& c, U, a)

        !This subroutine calculates the initial eikonal data S

        USE            omp_lib
        IMPLICIT       none

!--------------------------------------------------------!
!     Declare calling variables                          !
!--------------------------------------------------------!

        INTEGER*4               :: Nr, Nth, Nphi, Mr, Lmax, Lgrid
        
        REAL*8                  :: c
        REAL*8                  :: GLQWeights(Lgrid+1), GLQZeros(Lgrid+1)
        REAL*8                  :: rmax, rmin
        REAL*8                  :: U(Nth,Nphi)
        REAL*8                  :: rho(Nr)
        REAL*8                  :: theta(Nth)
        REAL*8                  :: phi(Nphi)

        COMPLEX*16, INTENT(out) :: a(Mr+1,2,Lmax+1,Lmax+1)

!--------------------------------------------------------!
!     Declare Locals                                     !
!--------------------------------------------------------!

        INTEGER*4        i, j, k
        INTEGER*4        crow      !Row counter
        COMPLEX*16       S(Nr,Nth,Nphi)
        REAL*8           r(Nr)

!--------------------------------------------------------!
!      Main Subroutine                                   !
!--------------------------------------------------------!

        PRINT *, 'Reinitializing!'

        CALL GetRadialCoordinates(Nr,rmax,rmin,rho,r)

        !$OMP PARALLEL DO SHARED(S, U, c) PRIVATE(i, k)
        DO j = 1, Nth
            DO k = 1, Nphi
                DO i = 1, Nr

                    S(i,j,k) = DCMPLX(100.0D0*(1.0D0+TANH((r(i)-U(j,k))/c)),&
                            &0.0D0 )

              END DO
           END DO
        END DO
        !$OMP END PARALLEL DO

        CALL SpatialToSpectralTransform(Nr,Nth,Nphi,Mr,Lmax,Lgrid,&
                                       &GLQWeights,GLQZeros,&
                                       &rho,theta,phi,&
                                       &S,a )

        RETURN
      END SUBROUTINE ReinitializeData
