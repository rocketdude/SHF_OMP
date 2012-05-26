!--------------------------------------------------------!
!    Evaluate Matrix of Metric subroutine                !
!--------------------------------------------------------!

      SUBROUTINE EvaluateMatrixofMetric(&
                 &alpha, beta1, beta2, beta3,&
                 &g11, g22, g33,&
                 &g12, g13, g23,&
                 &gAB)
                
        !This subroutine builds the matrix g_{A,B} given
        !alpha, beta^i, g_{ij}

        USE            omp_lib
        IMPLICIT       none

!--------------------------------------------------------!
!     Declare calling variables                          !
!--------------------------------------------------------!

        REAL*8, INTENT(in)   :: alpha
        REAL*8, INTENT(in)   :: beta1, beta2, beta3
        REAL*8, INTENT(in)   :: g11, g22, g33
        REAL*8, INTENT(in)   :: g12, g13, g23
        REAL*8, INTENT(out)  :: gAB(4,4)

!--------------------------------------------------------!
!     Declare local variables                            !
!--------------------------------------------------------!

        REAL*8                  det         !determinant of the spatial metric
        REAL*8                  gammaij     

!--------------------------------------------------------!
!      Main Subroutine                                   !
!--------------------------------------------------------!

        det
        gAB(1,1) = -1.0D0 * (alpha**2) + beta1
        

        RETURN
      END SUBROUTINE EvaluateMatrixofMetric
