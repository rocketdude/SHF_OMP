!----------------------------------------------------------!
!     GetMetricComponent  subroutine                       !
!----------------------------------------------------------!

    SUBROUTINE GetMetricComponent(&
    &Nr, Nth, Nphi,&
    &bufsize,&
    &time,&
    &iter, nchunks,&
    &DATASETFLAG,&
    &Filename,&
    &rmaxX, rmaxY, rmaxZ,&
    &rminX, rminY, rminZ,&
    &rho, theta, phi,&
    &MetricData)

        !We read in the metric data from the HDF5 file for a single iteration,
        !then interpolate the metric data which is in
        !Cartesian coordinates to match the spherical coordinates by
        !using Lekien-Marsden tricubic interpolation.
        
        !NOTE: DATASETFLAG =    0 --> alp: ADMBASE
        !                       1 --> betax : ADMBASE
        !                       2 --> betay : ADMBASE
        !                       3 --> betaz : ADMBASE
        !                       4 --> gxx : ADMBASE
        !                       5 --> gyy : ADMBASE
        !                       6 --> gzz : ADMBASE
        !                       7 --> gxy : ADMBASE
        !                       8 --> gxz : ADMBASE
        !                       9 --> gyz : ADMBASE


        USE omp_lib
        USE HDF5
        USE DynMetricArray

        IMPLICIT  none

!-----------------------------------------------------------!
!       Declare calling variables                           !
!-----------------------------------------------------------!

        INTEGER*4                   :: Nr, Nth, Nphi
        INTEGER*4                   :: iter         !iteration number in the CarpetCode simulation (iter != it)
        INTEGER*4                   :: nchunks      !number of chunks
        INTEGER*4                   :: DATASETFLAG

        INTEGER(HSIZE_T)            :: bufsize(3)        

        REAL*8                      :: rmaxX, rmaxY, rmaxZ
        REAL*8                      :: rminX, rminY, rminZ
        REAL*8                      :: rho(Nr)
        REAL*8                      :: theta(Nth)
        REAL*8                      :: phi(Nphi)

        CHARACTER*32                :: Filename(nchunks)

        REAL*8, INTENT(out)         :: MetricData(Nr,Nth,Nphi)
        REAL*8, INTENT(out)         :: time

!-----------------------------------------------------------!
!       Declare local variables                             !
!-----------------------------------------------------------!

        REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: metric10
        REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: metric9
        REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: metric8
        REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: metric7
        INTEGER*4                                dims10(3)
        INTEGER*4                                dims9(3)
        INTEGER*4                                dims8(3)
        INTEGER*4                                dims7(3)

        CHARACTER*32                format_string2
        CHARACTER*32                format_string1
        CHARACTER*100               dataset10(nchunks)
        CHARACTER*100               dataset9(nchunks)
        CHARACTER*100               dataset8(nchunks)
        CHARACTER*100               dataset7(nchunks)

        INTEGER*4                   Ox10, Oy10, Oz10    !Index for the origin for rl=10
        INTEGER*4                   Ox9, Oy9, Oz9       !Index for the origin for rl=9
        INTEGER*4                   Ox8, Oy8, Oz8       !Index for the origin for rl=8
        INTEGER*4                   Ox7, Oy7, Oz7       !Index for the origin for rl=7

        REAL*8                      Xmin10, Ymin10, Zmin10
        REAL*8                      Xmax10, Ymax10, Zmax10
        REAL*8                      Xmin9, Ymin9, Zmin9
        REAL*8                      Xmax9, Ymax9, Zmax9
        REAL*8                      Xmin8, Ymin8, Zmin8
        REAL*8                      Xmax8, Ymax8, Zmax8
        REAL*8                      Xmin7, Ymin7, Zmin7
        REAL*8                      Xmax7, Ymax7, Zmax7
        REAL*8                      delta10(3)          !Spatial discretization for rl=10
        REAL*8                      delta9(3)           !Spatial discretization for rl=9
        REAL*8                      delta8(3)           !Spatial discretization for rl=8
        REAL*8                      delta7(3)           !Spatial discretization for rl=7
        REAL*8                      time10, time9, time8, time7


        !Used for tricubic interpolation
        LOGICAL                     Inside10
        LOGICAL                     Inside9
        LOGICAL                     Inside8
        LOGICAL                     Inside7
        REAL*8                      r(Nr)
        REAL*8                      rX, rY, rZ
        REAL*8                      rmax, rmin
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

        !Allocation status
        INTEGER*4                   error

!----------------------------------------------------------!
!      Main                                                !
!----------------------------------------------------------!

        !----------------------------------------------------------!
        !      Write out dataset name                              !
        !----------------------------------------------------------!
        DO cnum = 1, nchunks

            IF( (DATASETFLAG .GE. 1) .AND. (DATASETFLAG .LE. 3) ) THEN
                IF( (cnum .LE. 10) .AND. (iter .LT. 10) ) THEN
                    format_string2 = '(A18,I1,A14,I1)'
                    format_string1 = '(A18,I1,A13,I1)'
                ELSE IF( (cnum .LE. 10) .AND. (iter .GE. 10) .AND. (iter .LT. 100) ) THEN
                    format_string2 = '(A18,I2,A14,I1)'
                    format_string1 = '(A18,I2,A13,I1)'
                ELSE IF( (cnum .LE. 10) .AND. (iter .GE. 100) .AND. (iter .LT. 1000) ) THEN
                    format_string2 = '(A18,I3,A14,I1)'
                    format_string1 = '(A18,I3,A13,I1)'
                ELSE IF( (cnum .LE. 10) .AND. (iter .GE. 1000) .AND. (iter .LT. 10000) ) THEN
                    format_string2 = '(A18,I4,A14,I1)'
                    format_string1 = '(A18,I4,A13,I1)'
                ELSE IF( (cnum .LE. 10) .AND. (iter .GE. 10000) .AND. (iter .LT. 100000) ) THEN
                    format_string2 = '(A18,I5,A14,I1)'
                    format_string1 = '(A18,I5,A13,I1)'
                ELSE IF( (cnum .GT. 10) .AND. (cnum .LE. 100) .AND. (iter .LT. 10) ) THEN
                    format_string2 = '(A18,I1,A14,I2)'
                    format_string1 = '(A18,I1,A13,I2)'
                ELSE IF( (cnum .GT. 10) .AND. (cnum .LE. 100) .AND. &
                        & (iter .GE. 10) .AND. (iter .LT. 100) ) THEN
                    format_string2 = '(A18,I2,A14,I2)'
                    format_string1 = '(A18,I2,A13,I2)'
                ELSE IF( (cnum .GT. 10) .AND. (cnum .LE. 100) .AND. &
                        & (iter .GE. 100) .AND. (iter .LT. 1000) ) THEN
                    format_string2 = '(A18,I3,A14,I2)'
                    format_string1 = '(A18,I3,A13,I2)'
                ELSE IF( (cnum .GT. 10) .AND. (cnum .LE. 100) .AND. &
                        & (iter .GE. 1000) .AND. (iter .LT. 10000) ) THEN
                    format_string2 = '(A18,I4,A14,I2)'
                    format_string1 = '(A18,I4,A13,I2)'
                ELSE IF( (cnum .GT. 10) .AND. (cnum .LE. 100) .AND. &
                        & (iter .GE. 10000) .AND. (iter .LT. 100000) ) THEN
                    format_string2 = '(A18,I5,A14,I2)'
                    format_string1 = '(A18,I5,A13,I2)'
                ELSE IF( (cnum .GT. 100) .AND. (cnum .LE. 1000) .AND. (iter .LT. 10) ) THEN
                    format_string2 = '(A18,I1,A14,I3)'
                    format_string1 = '(A18,I1,A13,I3)'
                ELSE IF( (cnum .GT. 100) .AND. (cnum .LE. 1000) .AND. &
                        & (iter .GE. 10) .AND. (iter .LT. 100) ) THEN
                    format_string2 = '(A18,I2,A14,I3)'
                    format_string1 = '(A18,I2,A13,I3)'
                ELSE IF( (cnum .GT. 100) .AND. (cnum .LE. 1000) .AND. &
                        & (iter .GE. 100) .AND. (iter .LT. 1000) ) THEN
                    format_string2 = '(A18,I3,A14,I3)'
                    format_string1 = '(A18,I3,A13,I3)'
                ELSE IF( (cnum .GT. 100) .AND. (cnum .LE. 1000) .AND. &
                        & (iter .GE. 1000) .AND. (iter .LT. 10000) ) THEN
                    format_string2 = '(A18,I4,A14,I3)'
                    format_string1 = '(A18,I4,A13,I3)'
                ELSE IF( (cnum .GT. 100) .AND. (cnum .LE. 1000) .AND. &
                        & (iter .GE. 10000) .AND. (iter .LT. 100000) ) THEN
                    format_string2 = '(A18,I5,A14,I3)'
                    format_string1 = '(A18,I5,A13,I3)'
                END IF

                IF( DATASETFLAG .EQ. 1 ) THEN
                    WRITE(dataset9(cnum), format_string1) 'ADMBASE::betax it=',iter,' tl=0 rl=9 c=',(cnum-1)
                    WRITE(dataset8(cnum), format_string1) 'ADMBASE::betax it=',iter,' tl=0 rl=8 c=',(cnum-1)
                    WRITE(dataset7(cnum), format_string1) 'ADMBASE::betax it=',iter,' tl=0 rl=7 c=',(cnum-1)
                    WRITE(dataset10(cnum), format_string2) 'ADMBASE::betax it=',iter,' tl=0 rl=10 c=',(cnum-1)
                ELSE IF( DATASETFLAG .EQ. 2 ) THEN
                    WRITE(dataset9(cnum), format_string1) 'ADMBASE::betay it=',iter,' tl=0 rl=9 c=',(cnum-1)
                    WRITE(dataset8(cnum), format_string1) 'ADMBASE::betay it=',iter,' tl=0 rl=8 c=',(cnum-1)
                    WRITE(dataset7(cnum), format_string1) 'ADMBASE::betay it=',iter,' tl=0 rl=7 c=',(cnum-1)
                    WRITE(dataset10(cnum), format_string2) 'ADMBASE::betay it=',iter,' tl=0 rl=10 c=',(cnum-1)
                ELSE IF( DATASETFLAG .EQ. 3 ) THEN
                    WRITE(dataset9(cnum), format_string1) 'ADMBASE::betaz it=',iter,' tl=0 rl=9 c=',(cnum-1)
                    WRITE(dataset8(cnum), format_string1) 'ADMBASE::betaz it=',iter,' tl=0 rl=8 c=',(cnum-1)
                    WRITE(dataset7(cnum), format_string1) 'ADMBASE::betaz it=',iter,' tl=0 rl=7 c=',(cnum-1)
                    WRITE(dataset10(cnum), format_string2) 'ADMBASE::betaz it=',iter,' tl=0 rl=10 c=',(cnum-1)
                END IF
            ELSE
                IF( (cnum .LE. 10) .AND. (iter .LT. 10) ) THEN
                    format_string2 = '(A16,I1,A14,I1)'
                    format_string1 = '(A16,I1,A13,I1)'
                ELSE IF( (cnum .LE. 10) .AND. (iter .GE. 10) .AND. (iter .LT. 100) ) THEN
                    format_string2 = '(A16,I2,A14,I1)'
                    format_string1 = '(A16,I2,A13,I1)'
                ELSE IF( (cnum .LE. 10) .AND. (iter .GE. 100) .AND. (iter .LT. 1000) ) THEN
                    format_string2 = '(A16,I3,A14,I1)'
                    format_string1 = '(A16,I3,A13,I1)'
                ELSE IF( (cnum .LE. 10) .AND. (iter .GE. 1000) .AND. (iter .LT. 10000) ) THEN
                    format_string2 = '(A16,I4,A14,I1)'
                    format_string1 = '(A16,I4,A13,I1)'
                ELSE IF( (cnum .LE. 10) .AND. (iter .GE. 10000) .AND. (iter .LT. 100000) ) THEN
                    format_string2 = '(A16,I5,A14,I1)'
                    format_string1 = '(A16,I5,A13,I1)'
                ELSE IF( (cnum .GT. 10) .AND. (cnum .LE. 100) .AND. (iter .LT. 10) ) THEN
                    format_string2 = '(A16,I1,A14,I2)'
                    format_string1 = '(A16,I1,A13,I2)'
                ELSE IF( (cnum .GT. 10) .AND. (cnum .LE. 100) .AND. &
                        & (iter .GE. 10) .AND. (iter .LT. 100) ) THEN
                    format_string2 = '(A16,I2,A14,I2)'
                    format_string1 = '(A16,I2,A13,I2)'
                ELSE IF( (cnum .GT. 10) .AND. (cnum .LE. 100) .AND. &
                        & (iter .GE. 100) .AND. (iter .LT. 1000) ) THEN
                    format_string2 = '(A16,I3,A14,I2)'
                    format_string1 = '(A16,I3,A13,I2)'
                ELSE IF( (cnum .GT. 10) .AND. (cnum .LE. 100) .AND. &
                        & (iter .GE. 1000) .AND. (iter .LT. 10000) ) THEN
                    format_string2 = '(A16,I4,A14,I2)'
                    format_string1 = '(A16,I4,A13,I2)'
                ELSE IF( (cnum .GT. 10) .AND. (cnum .LE. 100) .AND. &
                        & (iter .GE. 10000) .AND. (iter .LT. 100000) ) THEN
                    format_string2 = '(A16,I5,A14,I2)'
                    format_string1 = '(A16,I5,A13,I2)'
                ELSE IF( (cnum .GT. 100) .AND. (cnum .LE. 1000) .AND. (iter .LT. 10) ) THEN
                    format_string2 = '(A16,I1,A14,I3)'
                    format_string1 = '(A16,I1,A13,I3)'
                ELSE IF( (cnum .GT. 100) .AND. (cnum .LE. 1000) .AND. &
                        & (iter .GE. 10) .AND. (iter .LT. 100) ) THEN
                    format_string2 = '(A16,I2,A14,I3)'
                    format_string1 = '(A16,I2,A13,I3)'
                ELSE IF( (cnum .GT. 100) .AND. (cnum .LE. 1000) .AND. &
                        & (iter .GE. 100) .AND. (iter .LT. 1000) ) THEN
                    format_string2 = '(A16,I3,A14,I3)'
                    format_string1 = '(A16,I3,A13,I3)'
                ELSE IF( (cnum .GT. 100) .AND. (cnum .LE. 1000) .AND. &
                        & (iter .GE. 1000) .AND. (iter .LT. 10000) ) THEN
                    format_string2 = '(A16,I4,A14,I3)'
                    format_string1 = '(A16,I4,A13,I3)'
                ELSE IF( (cnum .GT. 100) .AND. (cnum .LE. 1000) .AND. &
                        & (iter .GE. 10000) .AND. (iter .LT. 100000) ) THEN
                    format_string2 = '(A16,I5,A14,I3)'
                    format_string1 = '(A16,I5,A13,I3)'
                END IF
                
                IF( DATASETFLAG .EQ. 0 ) THEN
                    WRITE(dataset9(cnum), format_string1) 'ADMBASE::alp it=',iter,' tl=0 rl=9 c=',(cnum-1)
                    WRITE(dataset8(cnum), format_string1) 'ADMBASE::alp it=',iter,' tl=0 rl=8 c=',(cnum-1)
                    WRITE(dataset7(cnum), format_string1) 'ADMBASE::alp it=',iter,' tl=0 rl=7 c=',(cnum-1)
                    WRITE(dataset10(cnum), format_string2) 'ADMBASE::alp it=',iter,' tl=0 rl=10 c=',(cnum-1)
                ELSE IF( DATASETFLAG .EQ. 4 ) THEN
                    WRITE(dataset9(cnum), format_string1) 'ADMBASE::gxx it=',iter,' tl=0 rl=9 c=',(cnum-1)
                    WRITE(dataset8(cnum), format_string1) 'ADMBASE::gxx it=',iter,' tl=0 rl=8 c=',(cnum-1)
                    WRITE(dataset7(cnum), format_string1) 'ADMBASE::gxx it=',iter,' tl=0 rl=7 c=',(cnum-1)
                    WRITE(dataset10(cnum), format_string2) 'ADMBASE::gxx it=',iter,' tl=0 rl=10 c=',(cnum-1)
                ELSE IF( DATASETFLAG .EQ. 5 ) THEN
                    WRITE(dataset9(cnum), format_string1) 'ADMBASE::gyy it=',iter,' tl=0 rl=9 c=',(cnum-1)
                    WRITE(dataset8(cnum), format_string1) 'ADMBASE::gyy it=',iter,' tl=0 rl=8 c=',(cnum-1)
                    WRITE(dataset7(cnum), format_string1) 'ADMBASE::gyy it=',iter,' tl=0 rl=7 c=',(cnum-1)
                    WRITE(dataset10(cnum), format_string2) 'ADMBASE::gyy it=',iter,' tl=0 rl=10 c=',(cnum-1)
                ELSE IF( DATASETFLAG .EQ. 6 ) THEN
                    WRITE(dataset9(cnum), format_string1) 'ADMBASE::gzz it=',iter,' tl=0 rl=9 c=',(cnum-1)
                    WRITE(dataset8(cnum), format_string1) 'ADMBASE::gzz it=',iter,' tl=0 rl=8 c=',(cnum-1)
                    WRITE(dataset7(cnum), format_string1) 'ADMBASE::gzz it=',iter,' tl=0 rl=7 c=',(cnum-1)
                    WRITE(dataset10(cnum), format_string2) 'ADMBASE::gzz it=',iter,' tl=0 rl=10 c=',(cnum-1)
                ELSE IF( DATASETFLAG .EQ. 7 ) THEN
                    WRITE(dataset9(cnum), format_string1) 'ADMBASE::gxy it=',iter,' tl=0 rl=9 c=',(cnum-1)
                    WRITE(dataset8(cnum), format_string1) 'ADMBASE::gxy it=',iter,' tl=0 rl=8 c=',(cnum-1)
                    WRITE(dataset7(cnum), format_string1) 'ADMBASE::gxy it=',iter,' tl=0 rl=7 c=',(cnum-1)
                    WRITE(dataset10(cnum), format_string2) 'ADMBASE::gxy it=',iter,' tl=0 rl=10 c=',(cnum-1)
                ELSE IF( DATASETFLAG .EQ. 8 ) THEN
                    WRITE(dataset9(cnum), format_string1) 'ADMBASE::gxz it=',iter,' tl=0 rl=9 c=',(cnum-1)
                    WRITE(dataset8(cnum), format_string1) 'ADMBASE::gxz it=',iter,' tl=0 rl=8 c=',(cnum-1)
                    WRITE(dataset7(cnum), format_string1) 'ADMBASE::gxz it=',iter,' tl=0 rl=7 c=',(cnum-1)
                    WRITE(dataset10(cnum), format_string2) 'ADMBASE::gxz it=',iter,' tl=0 rl=10 c=',(cnum-1)
                ELSE IF( DATASETFLAG .EQ. 9 ) THEN
                    WRITE(dataset9(cnum), format_string1) 'ADMBASE::gyz it=',iter,' tl=0 rl=9 c=',(cnum-1)
                    WRITE(dataset8(cnum), format_string1) 'ADMBASE::gyz it=',iter,' tl=0 rl=8 c=',(cnum-1)
                    WRITE(dataset7(cnum), format_string1) 'ADMBASE::gyz it=',iter,' tl=0 rl=7 c=',(cnum-1)
                    WRITE(dataset10(cnum), format_string2) 'ADMBASE::gyz it=',iter,' tl=0 rl=10 c=',(cnum-1)
                END IF

            END IF
        END DO

        !----------------------------------------------------------!
        !      Read rl=10 data                                     !
        !----------------------------------------------------------!

        CALL ReadHDF5MetricData(&
                &Filename, dataset10,&
                &DATASETFLAG,&
                &nchunks,&
                &time10,&
                &bufsize,&
                &Xmin10, Ymin10, Zmin10,&
                &Xmax10, Ymax10, Zmax10,&
                &delta10,&
                &Ox10, Oy10, Oz10)

        dims10 = metricdims
        ALLOCATE( metric10(dims10(1), dims10(2), dims10(3)), STAT=error )
        IF( error .NE. 0 ) STOP "*** Trouble allocating ***"
        
        metric10 = hdfmetric

        DEALLOCATE( hdfmetric, STAT=error )
        IF( error .NE. 0 ) STOP "*** Trouble deallocating ***"        

        !----------------------------------------------------------!
        !      Read rl=9 data                                      !
        !----------------------------------------------------------!

        CALL ReadHDF5MetricData(&
                &Filename, dataset9,&
                &DATASETFLAG,&
                &nchunks,&
                &time9,&
                &bufsize,&
                &Xmin9, Ymin9, Zmin9,&
                &Xmax9, Ymax9, Zmax9,&
                &delta9,&
                &Ox9, Oy9, Oz9)

        dims9 = metricdims
        ALLOCATE( metric9(dims9(1), dims9(2), dims9(3)), STAT=error )
        IF( error .NE. 0 ) STOP "*** Trouble allocating ***"

        metric9 = hdfmetric

        DEALLOCATE( hdfmetric, STAT=error )
        IF( error .NE. 0 ) STOP "*** Trouble deallocating ***"

        !----------------------------------------------------------!
        !      Read rl=8 data                                      !
        !----------------------------------------------------------!

        CALL ReadHDF5MetricData(&
                &Filename, dataset8,&
                &DATASETFLAG,&
                &nchunks,&
                &time8,&
                &bufsize,&
                &Xmin8, Ymin8, Zmin8,&
                &Xmax8, Ymax8, Zmax8,&
                &delta8,&
                &Ox8, Oy8, Oz8)

        dims8 = metricdims
        ALLOCATE( metric8(dims8(1), dims8(2), dims8(3)), STAT=error )
        IF( error .NE. 0 ) STOP "*** Trouble allocating ***"

        metric8 = hdfmetric

        DEALLOCATE( hdfmetric, STAT=error )
        IF( error .NE. 0 ) STOP "*** Trouble deallocating ***" 

        !----------------------------------------------------------!
        !      Read rl=7 data                                      !
        !----------------------------------------------------------!

        CALL ReadHDF5MetricData(&
                &Filename, dataset7,&
                &DATASETFLAG,&
                &nchunks,&
                &time7,&
                &bufsize,&
                &Xmin7, Ymin7, Zmin7,&
                &Xmax7, Ymax7, Zmax7,&
                &delta7,&
                &Ox7, Oy7, Oz7)

        dims7 = metricdims
        ALLOCATE( metric7(dims7(1), dims7(2), dims7(3)), STAT=error )
        IF( error .NE. 0 ) STOP "*** Trouble allocating ***"

        metric7 = hdfmetric

        DEALLOCATE( hdfmetric, STAT=error )
        IF( error .NE. 0 ) STOP "*** Trouble deallocating ***" 

!!!        IF( (time10 .NE. time9) .OR. (time10 .NE. time8) ) THEN
!!!            STOP "***ERROR*** Time discrepancy"
!!!        ELSE
!!!            time = time10
!!!        END IF

        time = time10
        
        !----------------------------------------------------------!
        !      Interpolate metric using                            !
        !      Lekien-Marsden tricubic interpolation routine       !
        !----------------------------------------------------------!

        IF( ALLOCATED(metric10) .AND. ALLOCATED(metric9) &
                &  .AND. ALLOCATED(metric8) .AND. ALLOCATED(metric7) ) THEN
        !$OMP PARALLEL DO &
        !$OMP &PRIVATE(i, k, x, y, z, Inside10, Inside9, Inside8, Inside7,&
        !$OMP & dx, dy, dz, &
        !$OMP & ni, nj, nk, x0, y0, z0, cube, fatxyz, &
        !$OMP & rmax, rmin, r)
        DO j = 1, Nth
            DO k = 1, Nphi
           
                CALL EvaluateRadialExtent(rmaxX,rmaxY,rmaxZ,&
                                        &theta(j),phi(k),rmax)
                CALL EvaluateRadialExtent(rminX,rminY,rminZ,&
                                        &theta(j),phi(k),rmin)

                r = 0.5D0*( (rmax-rmin) + (rmax-rmin)*rho )

                DO i = 1, Nr

                    x = r(i) * SIN( theta(j) ) * COS( phi(k) )
                    y = r(i) * SIN( theta(j) ) * SIN( phi(k) )
                    z = r(i) * COS( theta(j) ) 
                    
                    Inside10 =   (ABS(x) .LE. (Xmax10-2*delta10(1))) .AND.&
                                &(ABS(y) .LE. (Ymax10-2*delta10(2))) .AND.&
                                &(ABS(z) .LE. (Zmax10-2*delta10(3)))

                    Inside9 =    (ABS(x) .LE. (Xmax9-2*delta9(1))) .AND.&
                                &(ABS(y) .LE. (Ymax9-2*delta9(2))) .AND.&
                                &(ABS(z) .LE. (Zmax9-2*delta9(3))) 

                    Inside8 =    (ABS(x) .LE. (Xmax8-2*delta8(1))) .AND.&
                                &(ABS(y) .LE. (Ymax8-2*delta8(2))) .AND.&
                                &(ABS(z) .LE. (Zmax8-2*delta8(3))) 

                    Inside7 =    (ABS(x) .LE. (Xmax7-2*delta7(1))) .AND.&
                                &(ABS(y) .LE. (Ymax7-2*delta7(2))) .AND.&
                                &(ABS(z) .LE. (Zmax7-2*delta7(3))) 
    
                    IF( Inside10 .EQV. .TRUE. )THEN

                        dx = delta10(1)
                        dy = delta10(2)
                        dz = delta10(3)

                        ni = IDINT( x/dx ) + Ox10
                        nj = IDINT( y/dy ) + Oy10
                        nk = IDINT( z/dz ) + Oz10

                        IF( x .LT. 0.0D0 ) ni = ni-1
                        IF( y .LT. 0.0D0 ) nj = nj-1
                        IF( z .LT. 0.0D0 ) nk = nk-1

                        x0 = (ni-Ox10)*dx
                        y0 = (nj-Oy10)*dy
                        z0 = (nk-Oz10)*dz

                        cube = metric10( (ni-1):(ni+2), (nj-1):(nj+2), (nk-1):(nk+2) )

                    ELSEIF( (Inside10 .EQV. .FALSE.) .AND. (Inside9 .EQV. .TRUE.) ) THEN

                        dx = delta9(1)
                        dy = delta9(2)
                        dz = delta9(3)

                        ni = IDINT( x/dx ) + Ox9
                        nj = IDINT( y/dy ) + Oy9
                        nk = IDINT( z/dz ) + Oz9

                        IF( x .LT. 0.0D0 ) ni = ni-1
                        IF( y .LT. 0.0D0 ) nj = nj-1
                        IF( z .LT. 0.0D0 ) nk = nk-1

                        x0 = (ni-Ox9)*dx
                        y0 = (nj-Oy9)*dy
                        z0 = (nk-Oz9)*dz

                        cube = metric9( (ni-1):(ni+2), (nj-1):(nj+2), (nk-1):(nk+2) )

                    ELSEIF( (Inside9 .EQV. .FALSE.) .AND. (Inside8 .EQV. .TRUE.) ) THEN

                        dx = delta8(1)
                        dy = delta8(2)
                        dz = delta8(3)

                        ni = IDINT( x/dx ) + Ox8
                        nj = IDINT( y/dy ) + Oy8
                        nk = IDINT( z/dz ) + Oz8

                        IF( x .LT. 0.0D0 ) ni = ni-1
                        IF( y .LT. 0.0D0 ) nj = nj-1
                        IF( z .LT. 0.0D0 ) nk = nk-1

                        x0 = (ni-Ox8)*dx
                        y0 = (nj-Oy8)*dy
                        z0 = (nk-Oz8)*dz

                        cube = metric8( (ni-1):(ni+2), (nj-1):(nj+2), (nk-1):(nk+2) )

                    ELSEIF( (Inside8 .EQV. .FALSE.) .AND. (Inside7 .EQV. .TRUE.) ) THEN

                        dx = delta7(1)
                        dy = delta7(2)
                        dz = delta7(3)

                        ni = IDINT( x/dx ) + Ox7
                        nj = IDINT( y/dy ) + Oy7
                        nk = IDINT( z/dz ) + Oz7

                        IF( x .LT. 0.0D0 ) ni = ni-1
                        IF( y .LT. 0.0D0 ) nj = nj-1
                        IF( z .LT. 0.0D0 ) nk = nk-1

                        x0 = (ni-Ox7)*dx
                        y0 = (nj-Oy7)*dy
                        z0 = (nk-Oz7)*dz

                        cube = metric7( (ni-1):(ni+2), (nj-1):(nj+2), (nk-1):(nk+2) )
                                                                              
                    ELSE
                        PRINT *, 'ERROR: The spherical grid is larger than the grid for the metric'
                        STOP
                    END IF

                    CALL TricubicInterpolation(&
                            &cube,&
                            &dx, dy, dz,&
                            &x, y, z,&
                            &x0, y0, z0,&
                            &fatxyz)

                    MetricData(i,j,k) = fatxyz

                END DO
            END DO
        END DO
        !$OMP END PARALLEL DO

        ELSE
        PRINT *, 'ERROR in reading the metric data'
        END IF

        DEALLOCATE( metric10 )
        DEALLOCATE( metric9 )
        DEALLOCATE( metric8 )
        DEALLOCATE( metric7 )

        RETURN
      END SUBROUTINE GetMetricComponent
