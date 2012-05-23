!----------------------------------------------------------!
!     ReadMetric  subroutine                               !
!----------------------------------------------------------!

      SUBROUTINE ReadMetric(&
&M, Mr, NP,&
&h10, h9,&
&nx10, ny10, nz10,&
&nx9, ny9, nz9,&
&Ox10, Oy10, Oz10,&
&Ox9, Oy9, Oz9,&
&Xmax10, Ymax10, Zmax10,&
&Xmax9, Ymax9, Zmax9,&
&r, theta, phi,&
&LM,&
&MetricData)

        !We read in the metric data from the HDF5 file,
        !then interpolate the metric data which is in
        !Cartesian coordinates to match the spherical coordinates by
        !using Lekien-Marsden tricubic interpolation.

        USE       omp_lib
        IMPLICIT  none

!-----------------------------------------------------------!
!       Declare calling variables                           !
!-----------------------------------------------------------!

        INTEGER*4, INTENT(in)  :: M, Mr, NP
        INTEGER*4, INTENT(in)  :: nx10, ny10, nz10
        INTEGER*4, INTENT(in)  :: nx9, ny9, nz9
        INTEGER*4, INTENT(in)  :: Ox10, Oy10, Oz10
        INTEGER*4, INTENT(in)  :: Ox9, Oy9, Oz9

        INTEGER*4, INTENT(in)  :: LM(64, 64)

        REAL*8, INTENT(in)   :: h10, h9
        REAL*8, INTENT(in)   :: Xmax10, Ymax10, Zmax10
        REAL*8, INTENT(in)   :: Xmax9, Ymax9, Zmax9

        REAL*8, INTENT(in)   :: r(Mr+1)
        REAL*8, INTENT(in)   :: theta(2*M)
        REAL*8, INTENT(in)   :: phi(2*M)

        REAL*8, INTENT(out)  :: MetricData(4*NP)

!-----------------------------------------------------------!
!       Declare local variables                             !
!-----------------------------------------------------------!

        INTEGER*4                 i, j, k
        INTEGER*4                 ni, nj, nk
        INTEGER*4                 crow    !Row counter
        INTEGER*4                 RLFLAG

        CHARACTER*32              format_string

        REAL*8                    temp10(nx10, ny10, nz10)
        REAL*8                    temp9(nx9, ny9, nz9)
        REAL*8                    x, y, z
        REAL*8                    x0, y0, z0
        REAL*8                    dx, dy, dz
        REAL*8                    fatxyz

        REAL*8                    cube(4, 4, 4)

!----------------------------------------------------------!
!      Main                                                !
!----------------------------------------------------------!

        !CALL subroutine to read the metric
    
        !Interpolate metric using Lekien-Marsden tricubic interpolation routine
        DO i =  1, (Mr+1)
            DO j = 1, (2*M
                DO k = 1, (2*M)
                
                    crow = (i-1)*(2*M)*(2*M) + (j-1)*(2*M) + k
                
                    x = r(i) * SIN( theta(j) ) * COS( phi(k) )
                    y = r(i) * SIN( theta(j) ) * SIN( phi(k) )
                    z = r(i) * COS( theta(j) ) 
                    
                    IF( (ABS(x) .LE. (Xmax10-2*h10)) .AND.&
                        &(ABS(y) .LE. (Ymax10-2*h10)) .AND.&
                        &(ABS(z) .LE. (Zmax10-2*h10)) )THEN

                        RLFLAG = 10
                        ni = INT( x/h10 ) + Ox10
                        nj = INT( y/h10 ) + Oy10
                        nk = INT( z/h10 ) + Oz10

                        dx = h10
                        dy = h10
                        dz = h10

                        cube = temp10( (ni-1):(ni+2), (nj-1):(nj+2), (nk-1):(nk+2) )

                    ELSEIF( (ABS(x) .GT. Xmax10) .AND.&
                            &(ABS(y) .GT. Ymax10) .AND.&
                            &(ABS(z) .GT. Zmax10) .AND.&
                            &(ABS(x) .LE. (Xmax9-2*h9)) .AND.&
                            &(ABS(y) .LE. (Ymax9-2*h9)) .AND.&
                            &(ABS(z) .LE. (Zmax9-2*h9)) ) THEN

                        RLFLAG = 9
                        ni = INT( x/h9 ) + Ox9
                        nj = INT( y/h9 ) + Oy9
                        nk = INT( z/h9 ) + Oz9

                        dx = h9
                        dy = h9
                        dz = h9

                        cube = temp9( (ni-1):(ni+2), (nj-1):(nj+2), (nk-1):(nk+2) )

                    ELSE
                        PRINT *, 'ERROR: The spherical grid is larger than the grid for the metric'
                        STOP
                    END IF

                    x0 = ni*dx
                    y0 = nj*dy
                    z0 = nk*dz

                    CALL InterpolateMetricLM(&
                            &LM, cube,&
                            &dx, dy, dz,&
                            &x, y, z,&
                            &x0, y0, z0,&
                            &fatxyz)

                    MetricData(crow) = fatxyz

                END DO
            END DO
        END DO

        RETURN
      END SUBROUTINE ReadMetric
