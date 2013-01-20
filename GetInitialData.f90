!--------------------------------------------------------!
!    GetInitialData subroutine                           !
!--------------------------------------------------------!

      SUBROUTINE GetInitialData(& 
& Nr, Nth, Nphi,&
& Mr, Lmax, Lgrid,&
& GLQWeights, GLQZeros,&
& c,&
& X0, Y0, Z0,&
& rmaxX, rmaxY, rmaxZ,&
& rminX, rminY, rminZ,&
& rho, theta, phi,&
& a)

        !This subroutine calculates the initial eikonal data S

        USE            omp_lib
        IMPLICIT       none

!--------------------------------------------------------!
!     Declare calling variables                          !
!--------------------------------------------------------!

        INTEGER*4               :: Nr, Nth, Nphi
        INTEGER*4               :: Mr, Lmax, Lgrid
        
        REAL*8                  :: c
        REAL*8                  :: X0, Y0, Z0
        REAL*8                  :: rmaxX, rmaxY, rmaxZ
        REAL*8                  :: rminX, rminY, rminZ
        REAL*8                  :: GLQWeights(Lgrid+1), GLQZeros(Lgrid+1)
        REAL*8                  :: rho(Nr) 
        REAL*8                  :: theta(Nth)
        REAL*8                  :: phi(Nphi)

        COMPLEX*16, INTENT(out) :: a(Mr+1, 2, Lmax+1, Lmax+1)

!--------------------------------------------------------!
!     Declare Locals                                     !
!--------------------------------------------------------!

        INTEGER*4        i, j, k
        INTEGER*4        n
        COMPLEX*16       S(Nr, Nth, Nphi)
        REAL*8           r0
        REAL*8           rmax,rmin
        REAL*8           r(Nr)

!--------------------------------------------------------!
!      Main Subroutine                                   !
!--------------------------------------------------------!

        PRINT *, 'Evaluating initial S'

        !$OMP PARALLEL DO SHARED(S, r, c) &
        !$OMP &PRIVATE(k,i,r0,rmax,rmin,r)
        DO j = 1, Nth
            DO k = 1, Nphi
                
                CALL EvaluateRadialExtent(rmaxX,rmaxY,rmaxZ,&
                                         &theta(j),phi(k),rmax)
                CALL EvaluateRadialExtent(rminX,rminY,rminZ,&
                                         &theta(j),phi(k),rmin)

                r = 0.5D0*( (rmax-rmin) + (rmax-rmin)*rho )

                CALL EvaluateRadialExtent(X0,Y0,Z0,theta(j),phi(k),r0)

                DO i = 1, Nr

                   S(i,j,k) = DCMPLX( 100.0D0 *( 1 + TANH( ( r(i) - r0 )/c ) ),&
                                 & 0.0D0 )

              END DO
           END DO
        END DO
        !$OMP END PARALLEL DO

        CALL SpatialToSpectralTransform(Nr,Nth,Nphi,Mr,Lmax,Lgrid,&
                                       &GLQWeights,GLQZeros,&
                                       &rho,theta,phi,&
                                       &S,a )
        
        PRINT *, 'DONE!'

        RETURN
      END SUBROUTINE GetInitialData
