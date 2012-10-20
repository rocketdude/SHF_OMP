!--------------------------------------------------------!
!    GetInitialData subroutine                           !
!--------------------------------------------------------!

      SUBROUTINE GetInitialData(& 
& Nr, Nth, Nphi,&
& Mr, Lmax,&
& c,&
& R0,&
& r, rho, theta, phi,&
& a)

        !This subroutine calculates the initial eikonal data S

        USE            omp_lib
        IMPLICIT       none

!--------------------------------------------------------!
!     Declare calling variables                          !
!--------------------------------------------------------!

        INTEGER*4               :: Nr, Nth, Nphi
        INTEGER*4               :: Mr, Lmax
        
        REAL*8                  :: c
        REAL*8                  :: R0
        REAL*8                  :: r(Nr)
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


!--------------------------------------------------------!
!      Main Subroutine                                   !
!--------------------------------------------------------!

        PRINT *, 'Evaluating initial S'

        !$OMP PARALLEL DO SHARED(S, r, R0, c) PRIVATE(j, k)
        DO i = 1, Nr
           DO j = 1, Nth
              DO k = 1, Nphi

                 S(i,j,k) = CMPLX( 100.0D0 *( 1 + TANH( ( r(i) - R0 )/c ) ), &
                                 & 0.0D0 )

              END DO
           END DO
        END DO
        !$OMP END PARALLEL DO

        CALL SpatialToSpectralTransform(Nr,Nth,Nphi,Mr,Lmax,&
                                       &rho,theta,phi,&
                                       &S,a )
        
        PRINT *, 'DONE!'

        RETURN
      END SUBROUTINE GetInitialData
