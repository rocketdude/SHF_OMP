!--------------------------------------------------------!
!    Evaluate f subroutine                               !
!--------------------------------------------------------!

      SUBROUTINE EvaluateF(&
&n, l, ml,&
&rho, theta, phi,&
&f)

        !This subroutine calculates f_{nlm} at a specific collocation point

        USE             omp_lib
        IMPLICIT        none

!--------------------------------------------------------!
!     Declare calling variables                          !
!--------------------------------------------------------!

        INTEGER*4, INTENT(in)  :: n, l, ml

        REAL*8, INTENT(in)   :: rho
        REAL*8, INTENT(in)   :: theta
        REAL*8, INTENT(in)   :: phi
        
        COMPLEX*16, INTENT(out)  :: f

        COMPLEX*16      SphHarmonicY

!--------------------------------------------------------!
!     Declare Locals                                     !
!--------------------------------------------------------!

        COMPLEX*16      Tn
        COMPLEX*16      Ylm
        REAL*8          const

!--------------------------------------------------------!
!      Main Subroutine                                   !
!--------------------------------------------------------!

        IF( n .EQ. 0 ) THEN
           const = 0.5D0
        ELSE
           const = 1.0D0
        END IF

        Tn = CMPLX( COS( DBLE(n)*ACOS(rho) ), 0.0D0 )
        Ylm = SphHarmonicY(l, ml, theta, phi)

        f = CMPLX(const, 0.0D0)*Tn*Ylm

        RETURN
      END SUBROUTINE EvaluateF
