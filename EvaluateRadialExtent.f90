!--------------------------------------------------------!
!    Evaluate radial extent subroutine                   !
!--------------------------------------------------------!

      SUBROUTINE EvaluateRadialExtent(X,Y,Z,theta,phi,r) 

        USE             omp_lib
        IMPLICIT        none

!--------------------------------------------------------!
!     Declare calling variables                          !
!--------------------------------------------------------!

        REAL*8               :: X,Y,Z
        REAL*8               :: theta
        REAL*8               :: phi
        REAL*8, INTENT(out)  :: r

!--------------------------------------------------------!
!      Main Subroutine                                   !
!--------------------------------------------------------!

        rmax = SQRT( ( X*SIN(theta)*COS(phi) )**2 +&
                    &( Y*SIN(theta)*SIN(phi) )**2 +&
                    &( Z*COS(theta) )**2 )

        RETURN
      END SUBROUTINE EvaluateRadialExtent
