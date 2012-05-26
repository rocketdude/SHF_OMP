!----------------------------------------------------------!
!     TricubicInterpolation  subroutine                    !
!----------------------------------------------------------!

        SUBROUTINE TricubicInterpolation(&
&LM, cube,&
&dx, dy, dz,&
&x, y, z,&
&x0, y0, z0,&
&fatxyz)

        !Tricubic interpolation routine as specified by
        !Lekien and Marsden (Int. J. Numer. Meth. Eng 2005, 63).
        !This particular routine only works if the grid is equally spaced
        !in each direction with space dx, dy, dz.
        !You need to supply the values of to be interpolated on the cube
        !plus the values around the cube (+1 in each +/- direction)
        !in order to compute the partial derivatives.
        !So, instead of having the cube(2,2,2), we have cube(4,4,4)

        USE omp_lib
        IMPLICIT none

!----------------------------------------------------------!
!     Declare calling variables                            !
!----------------------------------------------------------!

        INTEGER*4, INTENT(in)   :: LM(64,64)

        REAL*8, INTENT(in)      :: cube(4,4,4)
        REAL*8, INTENT(in)      :: dx, dy, dz
        REAL*8, INTENT(in)      :: x, y, z
        REAL*8, INTENT(in)      :: x0, y0, z0

        REAL*8, INTENT(out)      :: fatxyz

!-----------------------------------------------------------!
!       Declare local variables                             !
!-----------------------------------------------------------!
        
        REAL*8                      alp(64)
        REAL*8                      b(64)

        REAL*8                      f(8)
        REAL*8                      dfdx(8)
        REAL*8                      dfdy(8)
        REAL*8                      dfdz(8)
        REAL*8                      d2fdxdy(8)
        REAL*8                      d2fdxdz(8)
        REAL*8                      d2fdydz(8)
        REAL*8                      d3fdxdydz(8)

        INTEGER*4                   i, j, k
        INTEGER*4                   ijk

!----------------------------------------------------------!
!      Main                                                !
!----------------------------------------------------------!

        !Setup the cube to be interpolated
        !This cube is of unit length in each direction so the values of
        !the derivatives are normalized.

        !Put f(i,j,k) -> f(1-8)
        f(1) = cube(2, 2, 2)
        f(2) = cube(3, 2, 2)
        f(3) = cube(2, 3, 2)
        f(4) = cube(3, 3, 2)

        f(5) = cube(2, 2, 3)
        f(6) = cube(3, 2, 3)
        f(7) = cube(2, 3, 3)
        f(8) = cube(3, 3, 3)

        !Compute the partial derivatives using finite difference
        
        !First derivative: central difference
        !dfdx(i,j,k) = (f(i+1, j, k)-f(i-1, j, k)) / 2
        dfdx(1) = (cube(2+1,2,2) - cube(2-1,2,2))/2.0D0
        dfdx(2) = (cube(3+1,2,2) - cube(3-1,2,2))/2.0D0
        dfdx(3) = (cube(2+1,3,2) - cube(2-1,3,2))/2.0D0
        dfdx(4) = (cube(3+1,3,2) - cube(3-1,3,2))/2.0D0

        dfdx(5) = (cube(2+1,2,3) - cube(2-1,2,3))/2.0D0
        dfdx(6) = (cube(3+1,2,3) - cube(3-1,2,3))/2.0D0
        dfdx(7) = (cube(2+1,3,3) - cube(2-1,3,3))/2.0D0
        dfdx(8) = (cube(3+1,3,3) - cube(3-1,3,3))/2.0D0

        !dfdy(i,j,k) = (f(i, j+1, k)-f(i, j-1, k)) / 2
        dfdy(1) = (cube(2,2+1,2) - cube(2,2-1,2))/2.0D0
        dfdy(2) = (cube(3,2+1,2) - cube(3,2-1,2))/2.0D0
        dfdy(3) = (cube(2,3+1,2) - cube(2,3-1,2))/2.0D0
        dfdy(4) = (cube(3,3+1,2) - cube(3,3-1,2))/2.0D0

        dfdy(5) = (cube(2,2+1,3) - cube(2,2-1,3))/2.0D0
        dfdy(6) = (cube(3,2+1,3) - cube(3,2-1,3))/2.0D0
        dfdy(7) = (cube(2,3+1,3) - cube(2,3-1,3))/2.0D0
        dfdy(8) = (cube(3,3+1,3) - cube(3,3-1,3))/2.0D0

        !dfdz(i,j,k) = (f(i, j, k+1)-f(i, j, k-1)) / 2
        dfdz(1) = (cube(2,2,2+1) - cube(2,2,2-1))/2.0D0
        dfdz(2) = (cube(3,2,2+1) - cube(3,2,2-1))/2.0D0
        dfdz(3) = (cube(2,3,2+1) - cube(2,3,2-1))/2.0D0
        dfdz(4) = (cube(3,3,2+1) - cube(3,3,2-1))/2.0D0

        dfdz(5) = (cube(2,2,3+1) - cube(2,2,3-1))/2.0D0
        dfdz(6) = (cube(3,2,3+1) - cube(3,2,3-1))/2.0D0
        dfdz(7) = (cube(2,3,3+1) - cube(2,3,3-1))/2.0D0
        dfdz(8) = (cube(3,3,3+1) - cube(3,3,3-1))/2.0D0

        !Second derivative: central difference
        !d2fdxdy(i,j,k) = (f(i+1,j+1,k)-f(i+1,j-1,k)-f(i-1,j+1,k)+f(i-1,j-1,k))/4
        d2fdxdy(1) = (cube(2+1,2+1,2)-cube(2+1,2-1,2)-cube(2-1,2+1,2)+cube(2-1,2-1,2))/4.0D0
        d2fdxdy(2) = (cube(3+1,2+1,2)-cube(3+1,2-1,2)-cube(3-1,2+1,2)+cube(3-1,2-1,2))/4.0D0
        d2fdxdy(3) = (cube(2+1,3+1,2)-cube(2+1,3-1,2)-cube(2-1,3+1,2)+cube(2-1,3-1,2))/4.0D0
        d2fdxdy(4) = (cube(3+1,3+1,2)-cube(3+1,3-1,2)-cube(3-1,3+1,2)+cube(3-1,3-1,2))/4.0D0

        d2fdxdy(5) = (cube(2+1,2+1,3)-cube(2+1,2-1,3)-cube(2-1,2+1,3)+cube(2-1,2-1,3))/4.0D0
        d2fdxdy(6) = (cube(3+1,2+1,3)-cube(3+1,2-1,3)-cube(3-1,2+1,3)+cube(3-1,2-1,3))/4.0D0
        d2fdxdy(7) = (cube(2+1,3+1,3)-cube(2+1,3-1,3)-cube(2-1,3+1,3)+cube(2-1,3-1,3))/4.0D0
        d2fdxdy(8) = (cube(3+1,3+1,3)-cube(3+1,3-1,3)-cube(3-1,3+1,3)+cube(3-1,3-1,3))/4.0D0

        !d2fdxdz(i,j,k) = (f(i+1,j,k+1)-f(i+1,j,k-1)-f(i-1,j,k+1)+f(i-1,j,k-1))/4
        d2fdxdz(1) = (cube(2+1,2,2+1)-cube(2+1,2,2-1)-cube(2-1,2,2+1)+cube(2-1,2,2-1))/4.0D0
        d2fdxdz(2) = (cube(3+1,2,2+1)-cube(3+1,2,2-1)-cube(3-1,2,2+1)+cube(3-1,2,2-1))/4.0D0
        d2fdxdz(3) = (cube(2+1,3,2+1)-cube(2+1,3,2-1)-cube(2-1,3,2+1)+cube(2-1,3,2-1))/4.0D0
        d2fdxdz(4) = (cube(3+1,3,2+1)-cube(3+1,3,2-1)-cube(3-1,3,2+1)+cube(3-1,3,2-1))/4.0D0

        d2fdxdz(5) = (cube(2+1,2,3+1)-cube(2+1,2,3-1)-cube(2-1,2,3+1)+cube(2-1,2,3-1))/4.0D0
        d2fdxdz(6) = (cube(3+1,2,3+1)-cube(3+1,2,3-1)-cube(3-1,2,3+1)+cube(3-1,2,3-1))/4.0D0
        d2fdxdz(7) = (cube(2+1,3,3+1)-cube(2+1,3,3-1)-cube(2-1,3,3+1)+cube(2-1,3,3-1))/4.0D0
        d2fdxdz(8) = (cube(3+1,3,3+1)-cube(3+1,3,3-1)-cube(3-1,3,3+1)+cube(3-1,3,3-1))/4.0D0

        !d2fdydz(i,j,k) = (f(i,j+1,k+1)-f(i,j+1,k-1)-f(i,j-1,k+1)+f(i,j-1,k-1))/4
        d2fdydz(1) = (cube(2,2+1,2+1)-cube(2,2+1,2-1)-cube(2,2-1,2+1)+cube(2,2-1,2-1))/4.0D0
        d2fdydz(2) = (cube(3,2+1,2+1)-cube(3,2+1,2-1)-cube(3,2-1,2+1)+cube(3,2-1,2-1))/4.0D0
        d2fdydz(3) = (cube(2,3+1,2+1)-cube(2,3+1,2-1)-cube(2,3-1,2+1)+cube(2,3-1,2-1))/4.0D0
        d2fdydz(4) = (cube(3,3+1,2+1)-cube(3,3+1,2-1)-cube(3,3-1,2+1)+cube(3,3-1,2-1))/4.0D0

        d2fdydz(5) = (cube(2,2+1,3+1)-cube(2,2+1,3-1)-cube(2,2-1,3+1)+cube(2,2-1,3-1))/4.0D0
        d2fdydz(6) = (cube(3,2+1,3+1)-cube(3,2+1,3-1)-cube(3,2-1,3+1)+cube(3,2-1,3-1))/4.0D0
        d2fdydz(7) = (cube(2,3+1,3+1)-cube(2,3+1,3-1)-cube(2,3-1,3+1)+cube(2,3-1,3-1))/4.0D0
        d2fdydz(8) = (cube(3,3+1,3+1)-cube(3,3+1,3-1)-cube(3,3-1,3+1)+cube(3,3-1,3-1))/4.0D0

        !Third derivative: central difference
        !d3fdxdydz(i,j,k) = (f(i+1,j+1,k+1) + f(i+1,j-1,k-1) + f(i-1,j+1,k-1) + f(i-1,j-1,k+1)
        !                   -f(i+1,j+1,k-1) - f(i+1,j-1,k+1) - f(i-1,j+1,k+1) - f(i-1,j-1,k-1))/8

        d3fdxdydz(1) = (cube(2+1,2+1,2+1) + cube(2+1,2-1,2-1) + cube(2-1,2+1,2-1) + cube(2-1,2-1,2+1) &
                    & - cube(2+1,2+1,2-1) - cube(2+1,2-1,2+1) - cube(2-1,2+1,2+1) - cube(2-1,2-1,2-1))/8.0D0
        d3fdxdydz(2) = (cube(3+1,2+1,2+1) + cube(3+1,2-1,2-1) + cube(3-1,2+1,2-1) + cube(3-1,2-1,2+1) &
                    & - cube(3+1,2+1,2-1) - cube(3+1,2-1,2+1) - cube(3-1,2+1,2+1) - cube(3-1,2-1,2-1))/8.0D0
        d3fdxdydz(3) = (cube(2+1,3+1,2+1) + cube(2+1,3-1,2-1) + cube(2-1,3+1,2-1) + cube(2-1,3-1,2+1) &
                    & - cube(2+1,3+1,2-1) - cube(2+1,3-1,2+1) - cube(2-1,3+1,2+1) - cube(2-1,3-1,2-1))/8.0D0
        d3fdxdydz(4) = (cube(3+1,3+1,2+1) + cube(3+1,3-1,2-1) + cube(3-1,3+1,2-1) + cube(3-1,3-1,2+1) &
                    & - cube(3+1,3+1,2-1) - cube(3+1,3-1,2+1) - cube(3-1,3+1,2+1) - cube(3-1,3-1,2-1))/8.0D0

        d3fdxdydz(5) = (cube(2+1,2+1,3+1) + cube(2+1,2-1,3-1) + cube(2-1,2+1,3-1) + cube(2-1,2-1,3+1) &
                    & - cube(2+1,2+1,3-1) - cube(2+1,2-1,3+1) - cube(2-1,2+1,3+1) - cube(2-1,2-1,3-1))/8.0D0
        d3fdxdydz(6) = (cube(3+1,2+1,3+1) + cube(3+1,2-1,3-1) + cube(3-1,2+1,3-1) + cube(3-1,2-1,3+1) &
                    & - cube(3+1,2+1,3-1) - cube(3+1,2-1,3+1) - cube(3-1,2+1,3+1) - cube(3-1,2-1,3-1))/8.0D0
        d3fdxdydz(7) = (cube(2+1,3+1,3+1) + cube(2+1,3-1,3-1) + cube(2-1,3+1,3-1) + cube(2-1,3-1,3+1) &
                    & - cube(2+1,3+1,3-1) - cube(2+1,3-1,3+1) - cube(2-1,3+1,3+1) - cube(2-1,3-1,3-1))/8.0D0
        d3fdxdydz(8) = (cube(3+1,3+1,3+1) + cube(3+1,3-1,3-1) + cube(3-1,3+1,3-1) + cube(3-1,3-1,3+1) &
                    & - cube(3+1,3+1,3-1) - cube(3+1,3-1,3+1) - cube(3-1,3+1,3+1) - cube(3-1,3-1,3-1))/8.0D0

        !Store the values of f and its derivatives in b
        DO i = 1, 8
            b(0+i) = f(i)
            b(8+i) = dfdx(i)
            b(16+i) = dfdy(i)
            b(24+i) = dfdz(i)
            b(32+i) = d2fdxdy(i)
            b(40+i) = d2fdxdz(i)
            b(48+i) = d2fdydz(i)
            b(56+i) = d2fdxdydz(i)
        END DO
        
        !Calculate the interpolation coffecients alp
        alp = MATMUL(LM, b)

        !Normalize the values of x, y, and z
        x = (x-x0)/dx
        y = (y-y0)/dy
        z = (z-z0)/dz

        !Calculate the value of x at the desired location
        fatxyz = 0.0D0
        DO i = 0, 3
            DO j = 0, 3
                DO k = 0, 3
                    ijk = 1 + i + 4*j * 16*k
                    fatxyz = fatxyz + alp(ijk) * (x**i) * (y**j) * (z**k)
                END DO
            END DO
        END DO

        RETURN
        END SUBROUTINE InterpolateMetricLM
