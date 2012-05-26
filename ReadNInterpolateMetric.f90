!----------------------------------------------------------!
!     Read & Interpolate Metric  subroutine                !
!----------------------------------------------------------!

    SUBROUTINE ReadNInterpolateMetric(&
    &M, Mr, NP,&
    &LM,&
    &bufsize,&
    &iter, nchunks,&
    &DATASETFLAG,&
    &Filename10, Filename9,&
    &r, theta, phi,&
    &MetricData)

        !We read in the metric data from the HDF5 file for a single iteration,
        !then interpolate the metric data which is in
        !Cartesian coordinates to match the spherical coordinates by
        !using Lekien-Marsden tricubic interpolation.
        
        !NOTE: DATASETFLAG =    0 --> alpha : KRANC2BSSNCHIMATTER
        !                       1 --> beta1 : KRANC2BSSNCHIMATTER
        !                       2 --> beta2 : KRANC2BSSNCHIMATTER
        !                       3 --> beta3 : KRANC2BSSNCHIMATTER
        !                       4 --> gxx : ADMBASE
        !                       5 --> gyy : ADMBASE
        !                       6 --> gzz : ADMBASE
        !                       7 --> gxy : ADMBASE
        !                       8 --> gxz : ADMBASE
        !                       9 --> gyz : ADMBASE


        USE omp_lib
        USE HDF5
        USE h5lt
        IMPLICIT  none

!-----------------------------------------------------------!
!       Declare calling variables                           !
!-----------------------------------------------------------!

        INTEGER*4, INTENT(in)       :: M, Mr, NP
        INTEGER*4, INTENT(in)       :: LM(64, 64)
        INTEGER*4, INTENT(in)       :: iter         !iteration number in the CarpetCode simulation (iter != it)
        INTEGER*4, INTENT(in)       :: nchunks      !number of chunks
        INTEGER*4, INTENT(in)       :: DATASETFLAG

        INTEGER(HSIZE_T), INTENT(in):: bufsize(3)        

        REAL*8, INTENT(in)          :: r(Mr+1)
        REAL*8, INTENT(in)          :: theta(2*M)
        REAL*8, INTENT(in)          :: phi(2*M)

        CHARACTER*32, INTENT(in)    :: Filename10
        CHARACTER*32, INTENT(in)    :: Filename9

        REAL*8, INTENT(out)         :: MetricData(4*NP)

!-----------------------------------------------------------!
!       Declare local variables                             !
!-----------------------------------------------------------!

        REAL*8, DIMENSION(:,:,:), POINTER ::    metric10
        REAL*8, DIMENSION(:,:,:), POINTER ::    metric9
        INTEGER*4                               dims10(3)
        INTEGER*4                               dims9(3)

        CHARACTER*32                format_string10
        CHARACTER*32                format_string9
        CHARACTER*32                dataset10(nchunks)
        CHARACTER*32                dataset9(nchunks)
        
        INTEGER*4                   Ox10, Oy10, Oz10    !Index for the origin for rl=10
        INTEGER*4                   Ox9, Oy9, Oz9       !Index for the origin for rl=9

        REAL*8                      Xmin10, Ymin10, Zmin10
        REAL*8                      Xmax10, Ymax10, Zmax10
        REAL*8                      Xmin9, Ymin9, Zmin9
        REAL*8                      Xmax9, Ymax9, Zmax9
        REAL*8                      delta10(3)          !Spatial discretization for rl=10
        REAL*8                      delta9(3)           !Spatial discretization for rl=9

        !Used for tricubic interpolation
        REAL*8                      x, y, z
        REAL*8                      x0, y0, z0
        REAL*8                      dx, dy, dz
        REAL*8                      fatxyz
        REAL*8                      cube(4, 4, 4)
        INTEGER*4                   ni, nj, nk

        !DO loop counters
        INTEGER*4                   cnum    !chunk number
        INTEGER*4                   i, j, k
        INTEGER*4                   crow    !Row counter


!----------------------------------------------------------!
!      Main                                                !
!----------------------------------------------------------!

        !----------------------------------------------------------!
        !      Write out dataset name                              !
        !----------------------------------------------------------!
        DO cnum = 1, nchunks

            IF( (DATASETFLAG .GE. 0) .AND. (DATASETFLAG .LE. 3) ) THEN
                IF( (cnum .LE. 10) .AND. (iter .LT. 10) ) THEN
                    format_string10 = '(A30,I1,A14,I1)'
                    format_string9 = '(A30,I1,A13,I1)'
                ELSE IF( (cnum .LE. 10) .AND. (iter .GE. 10) .AND. (iter .LT. 100) ) THEN
                    format_string10 = '(A30,I2,A14,I1)'
                    format_string9 = '(A30,I2,A13,I1)'
                ELSE IF( (cnum .LE. 10) .AND. (iter .GE. 100) .AND. (iter .LT. 1000) ) THEN
                    format_string10 = '(A30,I3,A14,I1)'
                    format_string9 = '(A30,I3,A13,I1)'
                ELSE IF( (cnum .LE. 10) .AND. (iter .GE. 100) .AND. (iter .LT. 10000) ) THEN
                    format_string10 = '(A30,I4,A14,I1)'
                    format_string9 = '(A30,I4,A13,I1)'
                ELSE IF( (cnum .GT. 10) .AND. (iter .LT. 10) ) THEN
                    format_string10 = '(A30,I1,A14,I2)'
                    format_string9 = '(A30,I1,A13,I2)'
                ELSE IF( (cnum .GT. 10) .AND. (iter .GE. 10) .AND. (iter .LT. 100) ) THEN
                    format_string10 = '(A30,I2,A14,I2)'
                    format_string9 = '(A30,I2,A13,I2)'
                ELSE IF( (cnum .GT. 10) .AND. (iter .GE. 100) .AND. (iter .LT. 1000) ) THEN
                    format_string10 = '(A30,I3,A14,I2)'
                    format_string9 = '(A30,I3,A13,I2)'
                ELSE IF( (cnum .GT. 10) .AND. (iter .GE. 100) .AND. (iter .LT. 10000) ) THEN
                    format_string10 = '(A30,I4,A14,I2)'
                    format_string9 = '(A30,I4,A13,I2)'
                END IF

                IF( DATASETFLAG .EQ. 0 ) THEN
                    WRITE(dataset9(cnum), format_string9) 'KRANC2BSSNCHIMATTER::alpha it=',iter,' tl=0 rl=9 c=',(cnum-1)
                    WRITE(dataset10(cnum), format_string10) 'KRANC2BSSNCHIMATTER::alpha it=',iter,' tl=0 rl=10 c=',(cnum-1)
                ELSE IF( DATASETFLAG .EQ. 1 ) THEN
                    WRITE(dataset9(cnum), format_string9) 'KRANC2BSSNCHIMATTER::beta1 it=',iter,' tl=0 rl=9 c=',(cnum-1)
                    WRITE(dataset10(cnum), format_string10) 'KRANC2BSSNCHIMATTER::beta1 it=',iter,' tl=0 rl=10 c=',(cnum-1)
                ELSE IF( DATASETFLAG .EQ. 2 ) THEN
                    WRITE(dataset9(cnum), format_string9) 'KRANC2BSSNCHIMATTER::beta2 it=',iter,' tl=0 rl=9 c=',(cnum-1)
                    WRITE(dataset10(cnum), format_string10) 'KRANC2BSSNCHIMATTER::beta2 it=',iter,' tl=0 rl=10 c=',(cnum-1)
                ELSE IF( DATASETFLAG .EQ. 3 ) THEN
                    WRITE(dataset9(cnum), format_string9) 'KRANC2BSSNCHIMATTER::beta3 it=',iter,' tl=0 rl=9 c=',(cnum-1)
                    WRITE(dataset10(cnum), format_string10) 'KRANC2BSSNCHIMATTER::beta3 it=',iter,' tl=0 rl=10 c=',(cnum-1)
                END IF
            ELSE
                IF( (cnum .LE. 10) .AND. (iter .LT. 10) ) THEN
                    format_string10 = '(A16,I1,A14,I1)'
                    format_string9 = '(A16,I1,A13,I1)'
                ELSE IF( (cnum .LE. 10) .AND. (iter .GE. 10) .AND. (iter .LT. 100) ) THEN
                    format_string10 = '(A16,I2,A14,I1)'
                    format_string9 = '(A16,I2,A13,I1)'
                ELSE IF( (cnum .LE. 10) .AND. (iter .GE. 100) .AND. (iter .LT. 1000) ) THEN
                    format_string10 = '(A16,I3,A14,I1)'
                    format_string9 = '(A16,I3,A13,I1)'
                ELSE IF( (cnum .LE. 10) .AND. (iter .GE. 100) .AND. (iter .LT. 10000) ) THEN
                    format_string10 = '(A16,I4,A14,I1)'
                    format_string9 = '(A16,I4,A13,I1)'
                ELSE IF( (cnum .GT. 10) .AND. (iter .LT. 10) ) THEN
                    format_string10 = '(A16,I1,A14,I2)'
                    format_string9 = '(A16,I1,A13,I2)'
                ELSE IF( (cnum .GT. 10) .AND. (iter .GE. 10) .AND. (iter .LT. 100) ) THEN
                    format_string10 = '(A16,I2,A14,I2)'
                    format_string9 = '(A16,I2,A13,I2)'
                ELSE IF( (cnum .GT. 10) .AND. (iter .GE. 100) .AND. (iter .LT. 1000) ) THEN
                    format_string10 = '(A16,I3,A14,I2)'
                    format_string9 = '(A16,I3,A13,I2)'
                ELSE IF( (cnum .GT. 10) .AND. (iter .GE. 100) .AND. (iter .LT. 10000) ) THEN
                    format_string10 = '(A16,I4,A14,I2)'
                    format_string9 = '(A16,I4,A13,I2)'
                END IF
                
                IF( DATASETFLAG .EQ. 4 ) THEN
                    WRITE(dataset9(cnum), format_string9) 'ADMBASE::gxx it=',iter,' tl=0 rl=9 c=',(cnum-1)
                    WRITE(dataset10(cnum), format_string10) 'ADMBASE::gxx it=',iter,' tl=0 rl=10 c=',(cnum-1)
                ELSE IF( DATASETFLAG .EQ. 5 ) THEN
                    WRITE(dataset9(cnum), format_string9) 'ADMBASE::gyy it=',iter,' tl=0 rl=9 c=',(cnum-1)
                    WRITE(dataset10(cnum), format_string10) 'ADMBASE::gyy it=',iter,' tl=0 rl=10 c=',(cnum-1)
                ELSE IF( DATASETFLAG .EQ. 6 ) THEN
                    WRITE(dataset9(cnum), format_string9) 'ADMBASE::gzz it=',iter,' tl=0 rl=9 c=',(cnum-1)
                    WRITE(dataset10(cnum), format_string10) 'ADMBASE::gzz it=',iter,' tl=0 rl=10 c=',(cnum-1)
                ELSE IF( DATASETFLAG .EQ. 7 ) THEN
                    WRITE(dataset9(cnum), format_string9) 'ADMBASE::gxy it=',iter,' tl=0 rl=9 c=',(cnum-1)
                    WRITE(dataset10(cnum), format_string10) 'ADMBASE::gxy it=',iter,' tl=0 rl=10 c=',(cnum-1)
                ELSE IF( DATASETFLAG .EQ. 8 ) THEN
                    WRITE(dataset9(cnum), format_string9) 'ADMBASE::gxz it=',iter,' tl=0 rl=9 c=',(cnum-1)
                    WRITE(dataset10(cnum), format_string10) 'ADMBASE::gxz it=',iter,' tl=0 rl=10 c=',(cnum-1)
                ELSE IF( DATASETFLAG .EQ. 9 ) THEN
                    WRITE(dataset9(cnum), format_string9) 'ADMBASE::gyz it=',iter,' tl=0 rl=9 c=',(cnum-1)
                    WRITE(dataset10(cnum), format_string10) 'ADMBASE::gyz it=',iter,' tl=0 rl=10 c=',(cnum-1)
                END IF

            END IF
        END DO

        !----------------------------------------------------------!
        !      Read rl=10 data                                     !
        !----------------------------------------------------------!

        CALL ReadHDF5Metric(&
                &Filename10, dataset10,&
                &nchunks,&
                &bufsize,&
                &Xmin10, Ymin10, Zmin10,&
                &Xmax10, Ymax10, Zmax10,&
                &delta10,&
                &Ox10, Oy10, Oz10,&
                &dims10,&
                &metric10)

        !----------------------------------------------------------!
        !      Read rl=9 data                                      !
        !----------------------------------------------------------!

        CALL ReadHDF5Metric(&
                &Filename9, dataset9,&
                &nchunks,&
                &bufsize,&
                &Xmin9, Ymin10, Zmin10,&
                &Xmax9, Ymax9, Zmax9,&
                &delta9,&
                &Ox9, Oy9, Oz9,&
                &dims9,&
                &metric9)

        !----------------------------------------------------------!
        !      Interpolate metric using                            !
        !      Lekien-Marsden tricubic interpolation routine       !
        !----------------------------------------------------------!

        DO i =  1, (Mr+1)
            DO j = 1, (2*M
                DO k = 1, (2*M)
                
                    crow = (i-1)*(2*M)*(2*M) + (j-1)*(2*M) + k
                
                    x = r(i) * SIN( theta(j) ) * COS( phi(k) )
                    y = r(i) * SIN( theta(j) ) * SIN( phi(k) )
                    z = r(i) * COS( theta(j) ) 
                    
                    IF( (ABS(x) .LE. (Xmax10-2*delta10(1))) .AND.&
                        &(ABS(y) .LE. (Ymax10-2*delta10(2))) .AND.&
                        &(ABS(z) .LE. (Zmax10-2*delta10(3))) )THEN

                        dx = delta10(1)
                        dy = delta10(2)
                        dz = delta10(3)

                        ni = IDINT( x/dx ) + Ox10
                        nj = IDINT( y/dy ) + Oy10
                        nk = IDINT( z/dz ) + Oz10

                        cube = metric10( (ni-1):(ni+2), (nj-1):(nj+2), (nk-1):(nk+2) )

                    ELSEIF( (ABS(x) .GT. Xmax10) .AND.&
                            &(ABS(y) .GT. Ymax10) .AND.&
                            &(ABS(z) .GT. Zmax10) .AND.&
                            &(ABS(x) .LE. (Xmax9-2*delta9(1))) .AND.&
                            &(ABS(y) .LE. (Ymax9-2*delta9(2))) .AND.&
                            &(ABS(z) .LE. (Zmax9-2*delta9(3))) ) THEN

                        dx = delta9(1)
                        dy = delta9(2)
                        dz = delta9(3)

                        ni = IDINT( x/dx ) + Ox9
                        nj = IDINT( y/dy ) + Oy9
                        nk = IDINT( z/dz ) + Oz9

                        cube = metric9( (ni-1):(ni+2), (nj-1):(nj+2), (nk-1):(nk+2) )

                    ELSE
                        PRINT *, 'ERROR: The spherical grid is larger than the grid for the metric'
                        STOP
                    END IF

                    x0 = ni*dx
                    y0 = nj*dy
                    z0 = nk*dz

                    CALL TricubicInterpolation(&
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
      END SUBROUTINE ReadNInterpolateMetric
