!---------------------------------------------------------!
!    FUNCTION Factorial                                   !
!---------------------------------------------------------!

  REAL*8 FUNCTION FactorialRatio(l, AbsM)

    !This function calculates the ratio of factorials:
    !(l-abs(ml))!/(l+abs(ml))!
    !found in the normalization of Spherical Harmonics

    IMPLICIT       none

    !---------------------------------------------------------!
    !     Declare variables passed by function                !
    !---------------------------------------------------------!

    INTEGER*4        l, AbsM

    !---------------------------------------------------------!
    !     Declare Locals                                      !
    !---------------------------------------------------------!	

    INTEGER*4        i   ! Counter for DO loops

    !---------------------------------------------------------!
    !     Calculate Factorial                                 !
    !---------------------------------------------------------!

    FactorialRatio = 1.0D0

    DO i = (l - AbsM +1), (l + AbsM)

       FactorialRatio = FactorialRatio / DBLE(i)

    END DO

    RETURN
  END FUNCTION FactorialRatio
