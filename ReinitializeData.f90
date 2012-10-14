!--------------------------------------------------------!
!    CalculateInitialData subroutine                     !
!--------------------------------------------------------!

      SUBROUTINE ReinitializeData(& 
& Nr, Nth, Nphi, Mr, Mlm,&
& r, rho, theta, phi,&
& c, U, a)

        !This subroutine calculates the initial eikonal data S

        USE            omp_lib
        IMPLICIT       none

!--------------------------------------------------------!
!     Declare calling variables                          !
!--------------------------------------------------------!

        INTEGER*4               :: Nr, Nth, Nphi, Mr, Mlm
        
        REAL*8                  :: c
        REAL*8                  :: U(Nth,Nphi)
        REAL*8                  :: r(Nr), rho(Nr)
        REAL*8                  :: theta(Nth)
        REAL*8                  :: phi(Nphi)

        COMPLEX*16, INTENT(out) :: a(Mr+1,Mlm)

!--------------------------------------------------------!
!     Declare Locals                                     !
!--------------------------------------------------------!

        INTEGER*4        i, j, k
        INTEGER*4        crow      !Row counter
        COMPLEX*16       S(Nr,Nth,Nphi)


!--------------------------------------------------------!
!      Main Subroutine                                   !
!--------------------------------------------------------!

        !$OMP PARALLEL DO SHARED(S, r, U, c) PRIVATE(j, k)
        DO i = 1, Nr
           DO j = 1, Nth
              DO k = 1, Nphi

                 S(i,j,k) = CMPLX( 100.0D0 *( 1 + TANH(( r(i) - U(j,k))/c ) ),&
                            &0.0D0 )

              END DO
           END DO
        END DO
        !$OMP END PARALLEL DO

        CALL SpatialToSpectralTransform(Nr,Nth,Nphi,Mr,Mlm,rho,theta,phi,S,a)

        RETURN
      END SUBROUTINE ReinitializeData
