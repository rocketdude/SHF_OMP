!--------------------------------------------------------!
!    Evaluate g subroutine                               !
!--------------------------------------------------------!

      SUBROUTINE EvaluateG(&
&n, l, ml,&
&rho, theta, phi,&
&M, Mr, g)

        !This subroutine calculates f_{nlm} at a specific collocation point

        USE             omp_lib
        IMPLICIT        none


!--------------------------------------------------------!
!     Declare calling variables                          !
!--------------------------------------------------------!

        INTEGER*4, INTENT(in)  :: n, l, ml
        INTEGER*4, INTENT(in)  :: M, Mr

        REAL*8, INTENT(in)   :: rho
        REAL*8, INTENT(in)   :: theta
        REAL*8, INTENT(in)   :: phi
        
        COMPLEX*16, INTENT(out)  :: g

        COMPLEX*16      SphHarmonicY

!--------------------------------------------------------!
!     Declare Locals                                     !
!--------------------------------------------------------!

        COMPLEX*16      Tn
        COMPLEX*16      Ylm

        REAL*8          const
        REAL*8          qsum
        REAL*8          wi, wj
        REAL*8          PI

        INTEGER*4         q, cq

!--------------------------------------------------------!
!      Main Subroutine                                   !
!--------------------------------------------------------!

        PI = 3.141592653589793238462643383279502884197D0

        const = SQRT(4.0D0*PI)*SQRT(2.0D0/PI)/DBLE(M)
        
        !Build wi (weight in the Chebyshev-Gauss quadrature)
        wi = PI/DBLE(Mr+1)

        !Build wj (weight in the Driscoll-Healy quadrature)
        qsum = 0
        DO q = 0, (M-1)
           cq = 2*q +1
           qsum = qsum + SIN( DBLE(cq)*theta )/DBLE(cq)
        END DO

        wj = 2.0D0*SQRT(2.0D0)/DBLE(2*M) * SIN(theta) * qsum

        Tn = CMPLX(COS( n*ACOS(rho) ), 0.0D0 )
        Ylm = SphHarmonicY(l, ml, theta, phi)

        g = CMPLX(const*wi*wj, 0.0D0)*CONJG(Tn*Ylm)

        RETURN
      END SUBROUTINE EvaluateG
