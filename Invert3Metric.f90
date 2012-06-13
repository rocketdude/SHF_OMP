!--------------------------------------------------------!
!    Invert 3-Metric subroutine                          !
!--------------------------------------------------------!

      SUBROUTINE Invert3Metric(g11, g22, g33,&
                              &g12, g13, g23,&
                              &u11, u22, u33,&
                              &u12, u13, u23,&
                              &error)
                
        !This subroutine inverts the spatial metric g_ij -> g^ij

        USE            omp_lib
        IMPLICIT       none

!--------------------------------------------------------!
!     Declare calling variables                          !
!--------------------------------------------------------!

        REAL*8                  :: g11, g22, g33 !g_ij
        REAL*8                  :: g12, g13, g23 !g_ij
        REAL*8, INTENT(out)     :: u11, u22, u33 !g^ij
        REAL*8, INTENT(out)     :: u12, u13, u23 !g^ij
        INTEGER*4, INTENT(out)  :: error         !error=0 if successful, =-1 if not

!--------------------------------------------------------!
!     Declare local variables                            !
!--------------------------------------------------------!

        REAL*8                  det

!--------------------------------------------------------!
!      Main Subroutine                                   !
!--------------------------------------------------------!

        error = 0

        !Calculate determinant
        det = g11*g22*g33 - g11*g23*g23&
            &-g12*g12*g33 + g12*g13*g23&
            &+g13*g12*g23 - g13*g13*g22
        
        IF( det .EQ. 0.0D0 ) error = -1
        
        u11 = (g22*g33 - g23*g23)/det
        u22 = (g11*g33 - g13*g13)/det
        u33 = (g11*g22 - g12*g12)/det
        u12 = (g13*g23 - g12*g33)/det
        u13 = (g12*g23 - g13*g22)/det
        u23 = (g12*g13 - g23*g11)/det

        RETURN
      END SUBROUTINE Invert3Metric
