!----------------------------------------------------------!
!     Evaluate Matrix of Metric  subroutine                !
!----------------------------------------------------------!

    SUBROUTINE EvaluateMatrixofMetric(&
               &alpha, beta1, beta2, beta3,&
               &g11, g22, g33,&
               &g12, g13, g23,&
               &gAB)

    !Build g with indices down
    !g's have indices down and beta's have indices up

    USE omp_lib
    IMPLICIT NONE

!-----------------------------------------------------------!
!       Declare calling variables                           !
!-----------------------------------------------------------!

    REAL*8, INTENT(in)                  :: alpha
    REAL*8, INTENT(in)                  :: beta1, beta2, beta3
    REAL*8, INTENT(in)                  :: g11, g22, g33
    REAL*8, INTENT(in)                  :: g12, g13, g23
    REAL*8, INTENT(out)                 :: gAB(4,4)

!-----------------------------------------------------------!
!       Main subroutine                                     !
!-----------------------------------------------------------!

    !spatial components
    gAB(2,2) = g11
    gAB(3,3) = g22
    gAB(4,4) = g33

    gAB(2,3) = g12
    gAB(3,2) = g12

    gAB(2,4) = g13
    gAB(4,2) = g13

    gAB(3,4) = g23
    gAB(4,3) = g23

    !beta_1 = beta^i * gamma_{i1}
    gAB(1,2) = beta1 * g11 + beta2 * g12 + beta3 * g13
    gAB(2,1) = gAB(1,2)

    !beta_2 = beta^i * gamma_{i2}
    gAB(1,3) = beta1 * g12 + beta2 * g22 + beta3 * g23
    gAB(3,1) = gAB(1,3)

    !beta_3 = beta^i * gamma_{i3}
    gAB(1,4) = beta1 * g13 + beta2 * g23 + beta3 * g33
    gAB(4,1) = gAB(1,4)

    !alpha
    gAB(1,1) = beta1*gAB(1,2) + beta2*gAB(1,3) + beta3*gAB(1,4) - ( alpha**2 )

    RETURN
    END SUBROUTINE EvaluateMatrixofMetric
