!--------------------------------------------------------!
!    GetInitialData subroutine                           !
!--------------------------------------------------------!

      SUBROUTINE GetInitialData(& 
& M, Mr, NP,&
& c,&
& R0,&
& r, theta, phi,&
& B,&
& a)

        !This subroutine calculates the initial eikonal data S

        USE            omp_lib
        IMPLICIT       none

!--------------------------------------------------------!
!     Declare calling variables                          !
!--------------------------------------------------------!

        INTEGER*4, INTENT(in)  :: M, Mr, NP
        
        REAL*8, INTENT(in)   :: c
        REAL*8, INTENT(in)   :: R0
        REAL*8, INTENT(in)   :: r(Mr+1)
        REAL*8, INTENT(in)   :: theta(2*M)
        REAL*8, INTENT(in)   :: phi(2*M)

        COMPLEX*16, INTENT(in)  :: B(NP, 4*NP)
        COMPLEX*16, INTENT(out) :: a(NP)

!--------------------------------------------------------!
!     Declare Locals                                     !
!--------------------------------------------------------!

        INTEGER*4        i, j, k
        INTEGER*4        crow      !Row counter
        COMPLEX*16     S(4*NP)


!--------------------------------------------------------!
!      Main Subroutine                                   !
!--------------------------------------------------------!

        PRINT *, 'Evaluating initial S'

        !$OMP PARALLEL DO SHARED(S, r, R0, c) PRIVATE(j, k, crow)
        DO i = 1, Mr+1
           DO j = 1, 2*M
              DO k = 1, 2*M

                 crow = (i-1)*(2*M)*(2*M) + (j-1)*(2*M) + k
                 S(crow) = CMPLX( 100.0D0 *( 1 + TANH( ( r(i) - R0 )/c ) ), 0.0D0 )

              END DO
           END DO
        END DO
        !$OMP END PARALLEL DO

        a = MATMUL(B, S)
        PRINT *, 'DONE!'

        RETURN
      END SUBROUTINE GetInitialData
