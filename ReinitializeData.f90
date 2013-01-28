!--------------------------------------------------------!
!    CalculateInitialData subroutine                     !
!--------------------------------------------------------!

      SUBROUTINE ReinitializeData(& 
& Nr, Nth, Nphi, &
& Mr, Lmax, Lgrid,&
& GLQWeights, GLQZeros,&
& rmaxX, rmaxY, rmaxZ,&
& rminX, rminY, rminZ,&
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
        REAL*8                  :: rmaxX, rmaxY, rmaxZ
        REAL*8                  :: rminX, rminY, rminZ
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
        REAL*8           rmax, rmin, r(Nr)

!--------------------------------------------------------!
!      Main Subroutine                                   !
!--------------------------------------------------------!

        PRINT *, 'Reinitializing!'
        !$OMP PARALLEL DO SHARED(S, U, c) PRIVATE(i, k, rmax, rmin, r)
        DO j = 1, Nth
            DO k = 1, Nphi

                CALL EvaluateRadialExtent(rmaxX,rmaxY,rmaxZ,&
                                         &theta(j),phi(k),rmax)
                CALL EvaluateRadialExtent(rminX,rminY,rminZ,&
                                         &theta(j),phi(k),rmin)
                r = 0.5D0*( (rmax-rmin) + (rmax-rmin)*rho )

                DO i = 1, Nr

                    S(i,j,k) = DCMPLX(100.0D0*(1+TANH(( r(i) - U(j,k))/c ) ),&
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
