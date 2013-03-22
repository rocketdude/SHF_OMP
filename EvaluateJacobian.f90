!--------------------------------------------------------!
!    Evaluate Jacobian subroutine                        !
!--------------------------------------------------------!

      SUBROUTINE EvaluateJacobian(&
                    &rmax,rmin,&
                    &rho,theta,phi,&
                    &JMatrix)
                
        !This subroutine builds the jacobian matrix to transform
        !the metric values in Cartesian coordinates (t,x,y,z)
        !to spherical coordinates (t,rho,theta,phi) at a specific collocation pt.

        USE             omp_lib
        IMPLICIT        none

!--------------------------------------------------------!
!     Declare calling variables                          !
!--------------------------------------------------------!

        REAL*8               :: rmax, rmin
        REAL*8               :: rho
        REAL*8               :: theta
        REAL*8               :: phi
        REAL*8, INTENT(out)  :: JMatrix(4,4)

!--------------------------------------------------------!
!     Declare local variables                            !
!--------------------------------------------------------!

        REAL*8               :: R

!--------------------------------------------------------!
!      Main Subroutine                                   !
!--------------------------------------------------------!

        R = 0.5D0*( (rmax+rmin)+(rmax-rmin)*rho )

        JMatrix(1,1) = 1.0D0
        JMatrix(1,2) = 0.0D0
        JMatrix(1,3) = 0.0D0
        JMatrix(1,4) = 0.0D0

        JMatrix(2,1) = 0.0D0
        JMatrix(2,2) = SIN(theta)*COS(phi)
        JMatrix(2,3) = R*COS(theta)*COS(phi)
        JMatrix(2,4) = -1.0D0*R*SIN(theta)*SIN(phi)

        JMatrix(3,1) = 0.0D0
        JMatrix(3,2) = SIN(theta)*SIN(phi)
        JMatrix(3,3) = R*COS(theta)*SIN(phi)
        JMatrix(3,4) = R*SIN(theta)*COS(phi)

        JMatrix(4,1) = 0.0D0
        JMatrix(4,2) = COS(theta)
        JMatrix(4,3) = -1.0D0*R*SIN(theta)
        JMatrix(4,4) = 0.0D0

        RETURN
      END SUBROUTINE EvaluateJacobian
!========================================================!
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

        r = SQRT( ( X*SIN(theta)*COS(phi) )**2 +&
                 &( Y*SIN(theta)*SIN(phi) )**2 +&
                 &( Z*COS(theta) )**2 )

        RETURN
      END SUBROUTINE EvaluateRadialExtent
!========================================================!
!--------------------------------------------------------!
!    Get radial coordinates                              !
!--------------------------------------------------------!

      SUBROUTINE GetRadialCoordinates(Nr,rmax,rmin,rho,r) 

        USE             omp_lib
        IMPLICIT        none

!--------------------------------------------------------!
!     Declare calling variables                          !
!--------------------------------------------------------!
        
        INTEGER*4            :: Nr
        REAL*8               :: rmax,rmin
        REAL*8               :: rho(Nr)
        REAL*8, INTENT(out)  :: r(Nr)

!--------------------------------------------------------!
!      Main Subroutine                                   !
!--------------------------------------------------------!

        r(:) = 0.5D0*( (rmax+rmin)+(rmax-rmin)*rho(:) )

        RETURN
      END SUBROUTINE GetRadialCoordinates
!========================================================!
!--------------------------------------------------------!
!    Get canonical radial coordinates                    !
!--------------------------------------------------------!

      SUBROUTINE GetRho(Nr,rmax,rmin,r,rho)

        USE             omp_lib
        IMPLICIT        none

!--------------------------------------------------------!
!     Declare calling variables                          !
!--------------------------------------------------------!

        INTEGER*4            :: Nr
        REAL*8               :: rmax,rmin
        REAL*8               :: r(Nr)
        REAL*8, INTENT(out)  :: rho(Nr)

!--------------------------------------------------------!
!      Main Subroutine                                   !
!--------------------------------------------------------!

        rho(:) = (2.0D0*r(:) - (rmax+rmin))/(rmax-rmin)

        RETURN
      END SUBROUTINE GetRho
