!--------------------------------------------------------!
!    Evaluate Jacobian subroutine                        !
!--------------------------------------------------------!

      SUBROUTINE EvaluateJacobian(r,theta,phi,JMatrix)
                
        !This subroutine builds the jacobian matrix to transform
        !the metric values in Cartesian coordinates
        !to spherical coordinates at a specific collocation point

        USE             omp_lib
        IMPLICIT        none

!--------------------------------------------------------!
!     Declare calling variables                          !
!--------------------------------------------------------!

        REAL*8               :: r
        REAL*8               :: theta
        REAL*8               :: phi
        REAL*8, INTENT(out)  :: JMatrix(4,4)

!--------------------------------------------------------!
!      Main Subroutine                                   !
!--------------------------------------------------------!

        JMatrix(1,1) = 1.0D0
        JMatrix(1,2) = 0.0D0
        JMatrix(1,3) = 0.0D0
        JMatrix(1,4) = 0.0D0

        JMatrix(2,1) = 0.0D0
        JMatrix(2,2) = SIN( theta ) * COS( phi )
        JMatrix(2,3) = r * COS( theta ) * COS( phi )
        JMatrix(2,4) = -1.0D0 * r * SIN( theta ) * SIN( phi )

        JMatrix(3,1) = 0.0D0
        JMatrix(3,2) = SIN( theta ) * SIN( phi )
        JMatrix(3,3) = r * COS( theta ) * SIN( phi )
        JMatrix(3,4) = r * SIN( theta ) * COS( phi )

        JMatrix(4,1) = 0.0D0
        JMatrix(4,2) = COS( theta )
        JMatrix(4,3) = -1.0D0 * r * SIN( theta )
        JMatrix(4,4) = 0.0D0

        RETURN
      END SUBROUTINE EvaluateJacobian
