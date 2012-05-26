!--------------------------------------------------------!
!    Evaluate Transformation Matrix subroutine           !
!--------------------------------------------------------!

      SUBROUTINE EvaluateTransformationMatrix(theta,phi,TMatrix)
                
        !This subroutine builds the transformation matrix from Cartesian coordinates
        !to spherical coordinates at a specific collocation point

        USE            omp_lib
        IMPLICIT       none

!--------------------------------------------------------!
!     Declare calling variables                          !
!--------------------------------------------------------!

        REAL*8, INTENT(in)   :: theta
        REAL*8, INTENT(in)   :: phi
        REAL*8, INTENT(out)  :: TMatrix(4,4)

!--------------------------------------------------------!
!      Main Subroutine                                   !
!--------------------------------------------------------!

        TMatrix(1,1) = 1.0D0
        TMatrix(1,2) = 0.0D0
        TMatrix(1,3) = 0.0D0
        TMatrix(1,4) = 0.0D0

        TMatrix(2,1) = 0.0D0
        TMatrix(2,2) = SIN( theta ) * COS( phi )
        TMatrix(2,3) = SIN( theta ) * SIN( phi )
        TMatrix(2,4) = COS( theta )

        TMatrix(3,1) = 0.0D0
        TMatrix(3,2) = COS( theta ) * COS( phi )
        TMatrix(3,3) = COS( theta ) * SIN( phi )
        TMatrix(3,4) = -1.0D0 * SIN( theta )

        TMatrix(4,1) = 0.0D0
        TMatrix(4,2) = -1.0D0 * SIN( phi )
        TMatrix(4,3) = COS( phi )
        TMatrix(4,4) = 0.0D0

        RETURN
      END SUBROUTINE EvaluateTransformationMatrix
