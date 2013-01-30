!--------------------------------------------------------!
!    Evaluate Jacobian subroutine                        !
!--------------------------------------------------------!

      SUBROUTINE EvaluateJacobian(&
                    &rmax,rmin,&
                    &drmaxdth,drmindth,drmaxdphi,drmindphi,&
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
        REAL*8               :: drmaxdth, drmindth
        REAL*8               :: drmaxdphi, drmindphi
        REAL*8               :: rho
        REAL*8               :: theta
        REAL*8               :: phi
        REAL*8, INTENT(out)  :: JMatrix(4,4)

!--------------------------------------------------------!
!     Declare local variables                            !
!--------------------------------------------------------!

        REAL*8               :: R,dR,dTHETA,dPHI

!--------------------------------------------------------!
!      Main Subroutine                                   !
!--------------------------------------------------------!

        R = 0.5D0*( (rmax+rmin)+(rmax-rmin)*rho )
        dR = 0.5D0*( rmax-rmin )
        dTHETA = 0.5D0*( (1.0D0+rho)*drmaxdth + (1.0D0-rho)*drmindth )
        dPHI = 0.5D0*( (1.0D0+rho)*drmaxdphi + (1.0D0-rho)*drmindphi )

        JMatrix(1,1) = 1.0D0
        JMatrix(1,2) = 0.0D0
        JMatrix(1,3) = 0.0D0
        JMatrix(1,4) = 0.0D0

        JMatrix(2,1) = 0.0D0
        JMatrix(2,2) = dR*SIN(theta)*COS(phi)
        JMatrix(2,3) = dTHETA*SIN(theta)*COS(phi)+R*COS(theta)*COS(phi)
        JMatrix(2,4) = dPHI*SIN(theta)*COS(phi)-R*SIN(theta)*SIN(phi)

        JMatrix(3,1) = 0.0D0
        JMatrix(3,2) = dR*SIN(theta)*SIN(phi)
        JMatrix(3,3) = dTHETA*SIN(theta)*SIN(phi)+R*COS(theta)*SIN(phi)
        JMatrix(3,4) = dPHI*SIN(theta)*SIN(phi)+R*SIN(theta)*COS(phi)

        JMatrix(4,1) = 0.0D0
        JMatrix(4,2) = dR*COS(theta)
        JMatrix(4,3) = dTHETA*COS(theta)-R*SIN(theta)
        JMatrix(4,4) = dPHI*COS(theta)

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
