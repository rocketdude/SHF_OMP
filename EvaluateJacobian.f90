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
!      Main Subroutine                                   !
!--------------------------------------------------------!

        JMatrix(1,1) = 1.0D0
        JMatrix(1,2) = 0.0D0
        JMatrix(1,3) = 0.0D0
        JMatrix(1,4) = 0.0D0

        JMatrix(2,1) = 0.0D0
        JMatrix(2,2) = 0.5D0 * (rmax - rmin) * SIN(theta) * COS(phi)
        JMatrix(2,3) = 0.5D0 * (rho + 1.0D0) * &
                        &( (drmaxdth - drmindth) * SIN(theta) * COS(phi) + &
                        &  (rmax - rmin) * COS(theta) * COS(phi) )
        JMatrix(2,4) = 0.5D0 * (rho + 1.0D0) * &
                        &( (drmaxdphi - drmindphi) * SIN(theta) * COS(phi) - &
                        &  (rmax - rmin) * SIN(theta) * SIN(phi) )

        JMatrix(3,1) = 0.0D0
        JMatrix(3,2) = 0.5D0 * (rmax - rmin) * SIN(theta) * SIN(phi)
        JMatrix(3,3) = 0.5D0 * (rho + 1.0D0) * &
                        &( (drmaxdth - drmindth) * SIN(theta) * SIN(phi) + &
                        &  (rmax - rmin) * COS(theta) * SIN(phi) )
        JMatrix(3,4) = 0.5D0 * (rho + 1.0D0) * &
                        &( (drmaxdphi - drmindphi) * SIN(theta) * SIN(phi) + &
                        &  (rmax - rmin) * SIN(theta) * COS(phi) )

        JMatrix(4,1) = 0.0D0
        JMatrix(4,2) = 0.5D0 * (rmax - rmin) * COS(theta)
        JMatrix(4,3) = 0.5D0 * (rho + 1.0D0) * &
                        &( (drmaxdth - drmindth) * COS(theta) - &
                        &  (rmax - rmin) * SIN(theta) )
        JMatrix(4,4) = 0.5D0 * (rho + 1.0D0) * &
                        &  (drmaxdphi - drmindphi) * COS(phi)

        RETURN
      END SUBROUTINE EvaluateJacobian
