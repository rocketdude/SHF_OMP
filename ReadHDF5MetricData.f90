!----------------------------------------------------------!
!     ReadHDF5MetricData  subroutine                       !
!----------------------------------------------------------!

    SUBROUTINE ReadHDF5MetricData(&
        &Filename, dataset,&
        &DATASETFLAG,&
        &CFLEN, CDLEN,&
        &nchunks,&
        &bufsize,&
        &Xmin, Ymin, Zmin,&
        &Xmax, Ymax, Zmax,&
        &delta,&
        &Ox, Oy, Oz)

        !This subroutine is to read the metric data supplied in HDF5 file.
        !Since the metric data is chunked into pieces, the subroutine
        !also recombines all the data into a single 3-dimensional array hdfmetric.
        !Our algorithm of recombining the metric data removes the outer buffer zones.
        !This algorithm assumes that Runge-Kutta 4 was used to evolve the data, such
        !that there are 4 substeps.

        !Uses Module: DynMetricArray

        USE omp_lib
        USE HDF5
        USE DynMetricArray

        IMPLICIT  none

!-----------------------------------------------------------!
!       Declare calling variables                           !
!-----------------------------------------------------------!

        INTEGER*4                           :: nchunks
        INTEGER*4                           :: CFLEN
        INTEGER*4                           :: CDLEN(nchunks)
        INTEGER*4                           :: DATASETFLAG
        INTEGER(HSIZE_T)                    :: bufsize(3)

        CHARACTER(LEN=*)                    :: Filename
        CHARACTER(LEN=*)                    :: dataset(nchunks)

        REAL*8, INTENT(out)                 :: Xmin, Ymin, Zmin
        REAL*8, INTENT(out)                 :: Xmax, Ymax, Zmax
        REAL*8, INTENT(out)                 :: delta(3)
        INTEGER*4, INTENT(out)              :: Ox, Oy, OZ

!-----------------------------------------------------------!
!       Declare local variables                             !
!-----------------------------------------------------------!

        !Used to read the HDF5 files
        INTEGER(HID_T)           :: file_id
        INTEGER(HID_T)           :: dset_id
        INTEGER(HID_T)           :: dspace_id
        INTEGER(HID_T)           :: attr_id      
        INTEGER(HSIZE_T)         :: nxyz(nchunks,3) !The dimensions of the dataset for each chunk
        INTEGER(HSIZE_T)         :: attr_dims(1)=3  !The dimensions of iorigin, origin, and delta are the same
        INTEGER*4                   hdferr  !Lets us know if there is error in reading the hdf5 files

        REAL*8                      buffer(nchunks,bufsize(1),bufsize(2),bufsize(3))
                                    !Buffer is the buffer for reading the dataset: buffer(chunk #, x, y, z)
                                    !Increase the buffer size if the dataset is larger, otherwise data will be truncated
        CHARACTER(LEN=CFLEN)            CFTemp
        CHARACTER(LEN=CDLEN(1))         CDTemp1     !This is a cheat, we're assuming there are only two different lengths of CDTemp
        CHARACTER(LEN=CDLEN(nchunks))   CDTemp2

        REAL*8, DIMENSION(:,:,:), ALLOCATABLE   :: combinedbuffer
        INTEGER*4                               :: combinedbufsize(3)

        !Used to chunk together the chunked pieces
        INTEGER*4                   ReflectWhichWay(3)  ! 0 --> no need to reflect
                                                        ! -1 --> reflect the bottom (negative) data to the top
                                                        ! +1 --> reflect the top (positive) data to the bottom
        LOGICAL                     PiSymmetry          ! .TRUE. if there's XY PiSymmetry
        INTEGER*4                   iorigin(nchunks,3)  !Bottom left index for each dataset
        INTEGER*4                   iTemp
        INTEGER(HSIZE_T)            minindex(3)
        INTEGER(HSIZE_T)            maxindex(3)
        INTEGER*4                   ghost(3)            !Extra grids due to ghost-zones
        INTEGER*4                   RKsubsteps
        INTEGER*4                   bufzone(3)
        REAL*8                      origin(nchunks,3)   !Bottom left values for each dataset
        REAL*8                      minpoint(3)
        REAL*8                      maxpoint(3)
        REAL*8                      fudge               !fudge factor

        !just placeholders in reading the HDF5 files--not really used
        INTEGER*4       ::          type_class
        INTEGER(SIZE_T) ::          type_size
        INTEGER(HSIZE_T)::          maxdims(3)

        !DO loop counters
        INTEGER*4                   i, j, k
        INTEGER*4                   ii, jj, kk
        INTEGER*4                   cnum                !chunk number

        !Allocation status
        INTEGER*4                   error

!----------------------------------------------------------!
!      Main                                                !
!----------------------------------------------------------!
 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !CHANGE THESE PARAMETERS WHEN NEEDED!
        fudge = 1.0D-7
        RKsubsteps = 4
        PiSymmetry = .TRUE.
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !Read the HDF5 files
        hdferr = 0

        !Initialize the Fortran interface
        CALL h5open_f(hdferr)

        !Open the file
        CFTemp = Filename
        CALL h5fopen_f(CFTemp, H5F_ACC_RDONLY_F, file_id, hdferr, H5P_DEFAULT_F)
        IF( hdferr .NE. 0 ) STOP "*** ERROR in opening file ***" 

        !Iterate through all the chunks
        IF( nchunks .LE. 10 ) THEN
        
            DO cnum = 1, nchunks
                CDTemp1 = dataset(cnum)
                
                !Open dataset
                CALL h5dopen_f(file_id, CDTemp1, dset_id, hdferr)
                IF( hdferr .NE. 0 ) STOP "*** ERROR in opening dataset ***" 

                !Get the dimensions
                CALL h5dget_space_f(dset_id, dspace_id, hdferr)
                IF( hdferr .NE. 0 ) STOP "*** ERROR in opening dataspace ***"
                CALL h5sget_simple_extent_dims_f(dspace_id, nxyz(cnum,:), maxdims, hdferr)
                IF( (nxyz(cnum,1) .GT. bufsize(1)) .OR. (nxyz(cnum,2) .GT. bufsize(2)) .OR. (nxyz(cnum,3) .GT. bufsize(3)) ) THEN
                    PRINT *, 'ERROR: buffer size is not big enough'
                    STOP
                END IF
                CALL h5sclose_f(dspace_id, hdferr)
                IF( hdferr .NE. 0 ) STOP "*** ERROR in closing dataspace ***"

                !Read dataset
                CALL h5dread_f(dset_id, H5T_IEEE_F64LE, buffer(cnum,1:nxyz(cnum,1),1:nxyz(cnum,2),1:nxyz(cnum,3)),&
                    &nxyz(cnum,:), hdferr, H5S_ALL_F, H5S_ALL_F, H5P_DEFAULT_F)
                IF( hdferr .NE. 0 ) STOP "*** ERROR in reading dataset ***"

                !Get the iorigin
                !CALL h5aopen_f(dset_id, 'iorigin', attr_id, hdferr)         !For HDF5 version 1.8.0 and beyond
                CALL h5aopen_name_f(dset_id, 'iorigin', attr_id, hdferr)    !For HDF5 version below 1.8.0
                IF( hdferr .NE. 0 ) STOP "*** ERROR in opening iorigin attribute***"
                CALL h5aread_f(attr_id, H5T_STD_I32LE, iorigin(cnum,:), attr_dims, hdferr)
                IF( hdferr .NE. 0 ) STOP "*** ERROR in getting iorigin ***"
                CALL h5aclose_f(attr_id, hdferr)
                IF( hdferr .NE. 0 ) STOP "*** ERROR in closing iorigin attribute***"

                !Get the values of origin
                !CALL h5aopen_f(dset_id, 'origin', attr_id, hdferr)         !For HDF5 version 1.8.0 and beyond
                CALL h5aopen_name_f(dset_id, 'origin', attr_id, hdferr)    !For HDF5 version below 1.8.0
                IF( hdferr .NE. 0 ) STOP "*** ERROR in opening origin attribute***"
                CALL h5aread_f(attr_id, H5T_IEEE_F64LE, origin(cnum,:), attr_dims, hdferr)
                IF( hdferr .NE. 0 ) STOP "*** ERROR in getting the values of origin ***"
                CALL h5aclose_f(attr_id, hdferr)
                IF( hdferr .NE. 0 ) STOP "*** ERROR in closing origin attribute***"
                
                IF( cnum .EQ. 1 ) THEN
                    !All the chunks should have the same spatial discretizations dx, dy, dz, 
                    !so just get these delta from a single dataset
                    !CALL h5aopen_f(dset_id, 'delta', attr_id, hdferr)         !For HDF5 version 1.8.0 and beyond
                    CALL h5aopen_name_f(dset_id, 'delta', attr_id, hdferr)    !For HDF5 version below 1.8.0
                    IF( hdferr .NE. 0 ) STOP "*** ERROR in opening delta attribute***"
                    CALL h5aread_f(attr_id, H5T_IEEE_F64LE, delta, attr_dims, hdferr)
                    IF( hdferr .NE. 0 ) STOP "*** ERROR in getting the values of dx, dy, dz ***"
                    CALL h5aclose_f(attr_id, hdferr)
                    IF( hdferr .NE. 0 ) STOP "*** ERROR in closing delta attribute***"
                END IF

                !Close dataset
                CALL h5dclose_f(dset_id, hdferr)
                IF( hdferr .NE. 0 ) STOP "*** ERROR in closing dataset ***" 

            END DO
        
        ELSE

            DO cnum = 1, 10
                CDTemp1 = dataset(cnum)
                
                !Open dataset
                CALL h5dopen_f(file_id, CDTemp1, dset_id, hdferr)
                IF( hdferr .NE. 0 ) STOP "*** ERROR in opening dataset ***" 

                !Get the dimensions
                CALL h5dget_space_f(dset_id, dspace_id, hdferr)
                IF( hdferr .NE. 0 ) STOP "*** ERROR in opening dataspace ***"
                CALL h5sget_simple_extent_dims_f(dspace_id, nxyz(cnum,:), maxdims, hdferr)
                IF( (nxyz(cnum,1) .GT. bufsize(1)) .OR. (nxyz(cnum,2) .GT. bufsize(2)) .OR. (nxyz(cnum,3) .GT. bufsize(3)) ) THEN
                    PRINT *, 'ERROR: buffer size is not big enough'
                    STOP
                END IF
                CALL h5sclose_f(dspace_id, hdferr)
                IF( hdferr .NE. 0 ) STOP "*** ERROR in closing dataspace ***"

                !Read dataset
                CALL h5dread_f(dset_id, H5T_IEEE_F64LE, buffer(cnum,1:nxyz(cnum,1),1:nxyz(cnum,2),1:nxyz(cnum,3)),&
                    &nxyz(cnum,:), hdferr, H5S_ALL_F, H5S_ALL_F, H5P_DEFAULT_F)
                IF( hdferr .NE. 0 ) STOP "*** ERROR in reading dataset ***"

                !Get the iorigin
                !CALL h5aopen_f(dset_id, 'iorigin', attr_id, hdferr)         !For HDF5 version 1.8.0 and beyond
                CALL h5aopen_name_f(dset_id, 'iorigin', attr_id, hdferr)    !For HDF5 version below 1.8.0
                IF( hdferr .NE. 0 ) STOP "*** ERROR in opening iorigin attribute***"
                CALL h5aread_f(attr_id, H5T_STD_I32LE, iorigin(cnum,:), attr_dims, hdferr)
                IF( hdferr .NE. 0 ) STOP "*** ERROR in getting iorigin ***"
                CALL h5aclose_f(attr_id, hdferr)
                IF( hdferr .NE. 0 ) STOP "*** ERROR in closing iorigin attribute***"

                !Get the values of origin
                !CALL h5aopen_f(dset_id, 'origin', attr_id, hdferr)         !For HDF5 version 1.8.0 and beyond
                CALL h5aopen_name_f(dset_id, 'origin', attr_id, hdferr)    !For HDF5 version below 1.8.0
                IF( hdferr .NE. 0 ) STOP "*** ERROR in opening origin attribute***"
                CALL h5aread_f(attr_id, H5T_IEEE_F64LE, origin(cnum,:), attr_dims, hdferr)
                IF( hdferr .NE. 0 ) STOP "*** ERROR in getting the values of origin ***"
                CALL h5aclose_f(attr_id, hdferr)
                IF( hdferr .NE. 0 ) STOP "*** ERROR in closing origin attribute***"

                IF( cnum .EQ. 1 ) THEN
                    !All the chunks should have the same spatial discretizations dx, dy, dz, 
                    !so just get these delta from a single dataset
                    !CALL h5aopen_f(dset_id, 'delta', attr_id, hdferr)         !For HDF5 version 1.8.0 and beyond
                    CALL h5aopen_name_f(dset_id, 'delta', attr_id, hdferr)    !For HDF5 version below 1.8.0
                    IF( hdferr .NE. 0 ) STOP "*** ERROR in opening delta attribute***"
                    CALL h5aread_f(attr_id, H5T_IEEE_F64LE, delta, attr_dims, hdferr)
                    IF( hdferr .NE. 0 ) STOP "*** ERROR in getting the values of dx, dy, dz ***"
                    CALL h5aclose_f(attr_id, hdferr)
                    IF( hdferr .NE. 0 ) STOP "*** ERROR in closing delta attribute***"
                END IF

                !Close dataset
                CALL h5dclose_f(dset_id, hdferr)
                IF( hdferr .NE. 0 ) STOP "*** ERROR in closing dataset ***" 

            END DO

            DO cnum = 11, nchunks
                CDTemp2 = dataset(cnum)

                !Open dataset
                CALL h5dopen_f(file_id, CDTemp2, dset_id, hdferr)
                IF( hdferr .NE. 0 ) STOP "*** ERROR in opening dataset ***" 

                !Get the dimensions
                CALL h5dget_space_f(dset_id, dspace_id, hdferr)
                IF( hdferr .NE. 0 ) STOP "*** ERROR in opening dataspace ***"
                CALL h5sget_simple_extent_dims_f(dspace_id, nxyz(cnum,:), maxdims, hdferr)
                IF( (nxyz(cnum,1) .GT. bufsize(1)) .OR. (nxyz(cnum,2) .GT. bufsize(2)) .OR. (nxyz(cnum,3) .GT. bufsize(3)) ) THEN
                    PRINT *, 'ERROR: buffer size is not big enough'
                    STOP
                END IF
                CALL h5sclose_f(dspace_id, hdferr)
                IF( hdferr .NE. 0 ) STOP "*** ERROR in closing dataspace ***"

                !Read dataset
                CALL h5dread_f(dset_id, H5T_IEEE_F64LE, buffer(cnum,1:nxyz(cnum,1),1:nxyz(cnum,2),1:nxyz(cnum,3)),&
                    &nxyz(cnum,:), hdferr, H5S_ALL_F, H5S_ALL_F, H5P_DEFAULT_F)
                IF( hdferr .NE. 0 ) STOP "*** ERROR in reading dataset ***"

                !Get the iorigin
                !CALL h5aopen_f(dset_id, 'iorigin', attr_id, hdferr)         !For HDF5 version 1.8.0 and beyond
                CALL h5aopen_name_f(dset_id, 'iorigin', attr_id, hdferr)    !For HDF5 version below 1.8.0
                IF( hdferr .NE. 0 ) STOP "*** ERROR in opening iorigin attribute***"
                CALL h5aread_f(attr_id, H5T_STD_I32LE, iorigin(cnum,:), attr_dims, hdferr)
                IF( hdferr .NE. 0 ) STOP "*** ERROR in getting iorigin ***"
                CALL h5aclose_f(attr_id, hdferr)
                IF( hdferr .NE. 0 ) STOP "*** ERROR in closing iorigin attribute***"

                !Get the values of origin
                !CALL h5aopen_f(dset_id, 'origin', attr_id, hdferr)         !For HDF5 version 1.8.0 and beyond
                CALL h5aopen_name_f(dset_id, 'origin', attr_id, hdferr)    !For HDF5 version below 1.8.0
                IF( hdferr .NE. 0 ) STOP "*** ERROR in opening origin attribute***"
                CALL h5aread_f(attr_id, H5T_IEEE_F64LE, origin(cnum,:), attr_dims, hdferr)
                IF( hdferr .NE. 0 ) STOP "*** ERROR in getting the values of origin ***"
                CALL h5aclose_f(attr_id, hdferr)
                IF( hdferr .NE. 0 ) STOP "*** ERROR in closing origin attribute***"

                !Close dataset
                CALL h5dclose_f(dset_id, hdferr)
                IF( hdferr .NE. 0 ) STOP "*** ERROR in closing dataset ***" 

            END DO

        END IF
        !Close the file
        CALL h5fclose_f(file_id, hdferr)
        IF( hdferr .NE. 0 ) STOP "*** ERROR in closing hdf file ***"
        CALL h5close_f( hdferr )

        !PROCESS THE DATA
        DO i = 1,3
            !Find the left-bottom corner (minimum index)
            minindex(i) = MINVAL(iorigin(:,i))

            !Find the top-right corner (maximum index)
            maxindex(i) = MAXVAL(iorigin(:,i)+nxyz(:,i))-1

            !Make minindex = (1,1,1) (left-bottom) and change the maxindex & iorigin accordingly
            iorigin(:,i) = iorigin(:,i) - minindex(i) +1
            maxindex(i) = maxindex(i) - minindex(i) +1
            minindex(i) = 1

            !Find the minimum x,y,and z
            minpoint(i) = MINVAL(origin(:,i))
            !Find the maximum x,y,and z
            maxpoint(i) = minpoint(i) + DBLE(maxindex(i)-minindex(i))*delta(i)

            !Find out if we either only have half the data or all the data, in each direction
            IF( ABS(maxpoint(i)) .LT. (ABS(minpoint(i))+fudge) .AND.&
              & ABS(maxpoint(i)) .GT. (ABS(minpoint(i))-fudge) ) THEN
                combinedbufsize(i) = maxindex(i)
                ReflectWhichWay(i) = 0
                
                !But now we need to calculate what the size of the ghost zone is
                ghost(i) = (iorigin(1,i) + nxyz(1,i) - iorigin(2,i))/2
                IF( nchunks .GE. 3 ) THEN
                    DO cnum = 3, nchunks
                        iTemp = (iorigin(1,i) + nxyz(1,i) - iorigin(cnum,i))/2
                        IF( (ghost(i) .EQ. 0) .AND. (iTemp .GT. 0) ) THEN
                            ghost(i) = iTemp
                        ELSE IF( (ghost(i) .NE. 0) .AND. (iTemp .GT. 0) .AND. (iTemp .LT. ghost(i)) ) THEN
                            ghost(i) = iTemp
                        END IF
                    END DO
                END IF

            ELSE IF( ABS(maxpoint(i)) .GT. (ABS(minpoint(i))+fudge) ) THEN
                IF( minpoint(i) .GT. 0 ) THEN
                    PRINT *, 'ERROR: in spatial index', i, ', the minimum point of the metric has to be less than 0'
                    STOP
                ELSE
                    ghost(i) = IDINT(fudge + ABS(minpoint(i))/delta(i) ) !Any data below the 0 plane would be ghost zone
                    maxindex(i) = maxindex(i) - ghost(i)    !remove the bottom ghost zones but keep minindex = 1
                    combinedbufsize(i) = 2*( maxindex(i)-minindex(i) ) + 1
                    ReflectWhichWay(i) = +1
                END IF
            ELSE IF( ABS(maxpoint(i)) .LT. (ABS(minpoint(i))-fudge) ) THEN
                IF( maxpoint(i) .LT. 0 ) THEN
                    PRINT *, 'ERROR: in spatial index', i, ', the maximum point of the metric has to be more than 0'
                    STOP
                ELSE
                    ghost(i) = IDINT(fudge + maxpoint(i)/delta(i))  !Any data above the 0 plane would be ghost zone
                    maxindex(i) = maxindex(i) - ghost(i)     !remove the top ghost zones and keep minindex = 1
                    combinedbufsize(i) = 2*( maxindex(i)-minindex(i) ) + 1 
                    ReflectWhichWay(i) = -1
                END IF
            END IF

        END DO

        !Allocate combinedbuffer to its appropriate size
        ALLOCATE( combinedbuffer(combinedbufsize(1), combinedbufsize(2), combinedbufsize(3)), STAT=error )
        IF( error .NE. 0 ) STOP "*** Trouble Allocating ***"

        !Combine the chunked data and save it in combinedbuffer
        !$OMP PARALLEL DO PRIVATE(i,j,k,ii,jj,kk)
        DO cnum = 1, nchunks
            DO i = 1, nxyz(cnum,1)
                DO j = 1, nxyz(cnum,2)
                    DO k = 1, nxyz(cnum,3)
                        ii = (i+iorigin(cnum,1)-1)
                        jj = (j+iorigin(cnum,2)-1)
                        kk = (k+iorigin(cnum,3)-1)
                        combinedbuffer(ii,jj,kk) = buffer(cnum,i,j,k)   !NOTE:combinedbuffer still includes ghost zones
                                                                        !     and buffer zones
                    END DO
                END DO
            END DO
        END DO
        !$OMP END PARALLEL DO

        !Perform the reflection accordingly
        !NOTE: There is no direct relationship between (i,j,k)/(ii,jj,kk) and the coordinates (x,y,z)
        !(i,j,k) and (ii,jj,kk) are really used as dummy place holders here
        
        !X-DIRECTION
        IF( ReflectWhichWay(1) .EQ. +1 ) THEN !DATA IS FOR THE TOP HEMISPHERE
            !The data is now located at the bottom hemisphere, move it to the top hemisphere, remove ghost zones
            DO i = minindex(1), maxindex(1)
                jj = combinedbufsize(1) -i+1
                kk = maxindex(1) + ghost(1) -i+1
                combinedbuffer(jj,:,:) = combinedbuffer(kk,:,:)
            END DO
            !Reflect the top hemisphere values about x = 0
            DO i = minindex(1), (maxindex(1)-1)
                jj = combinedbufsize(1) -i +1
                combinedbuffer(i,:,:) = combinedbuffer(jj,:,:)
            END DO
        ELSE IF( ReflectWhichWay(1) .EQ. -1 ) THEN !DATA IS FOR THE BOTTOM HEMISPHERE
            !Reflect the bottom hemisphere values about x = 0
            DO i = minindex(1), (maxindex(1)-1)
                jj = combinedbufsize(1) -i +1
                combinedbuffer(jj,:,:) = combinedbuffer(i,:,:)
            END DO
        END IF

        !Y-DIRECTION
        IF( ReflectWhichWay(2) .EQ. +1 ) THEN !DATA IS FOR THE TOP HEMISPHERE
           !The data is now located at the bottom hemisphere, move it to the top hemisphere, remove ghost zones
            DO i = minindex(2), maxindex(2)
                jj = combinedbufsize(2) -i+1
                kk = maxindex(2) + ghost(2) -i+1
                combinedbuffer(:,jj,:) = combinedbuffer(:,kk,:)
            END DO
            !Reflect the top hemisphere values about y = 0
            DO i = minindex(2), (maxindex(2)-1)
                jj = combinedbufsize(2) -i +1
                combinedbuffer(:,i,:) = combinedbuffer(:,jj,:)
            END DO
        ELSE IF( ReflectWhichWay(2) .EQ. -1 ) THEN !DATA IS FOR THE BOTTOM HEMISPHERE
            !Reflect the bottom hemisphere values about y = 0
            DO i = minindex(2), (maxindex(2)-1)
                jj = combinedbufsize(2) -i +1
                combinedbuffer(:,jj,:) = combinedbuffer(:,i,:)
            END DO
        END IF

        !Check for PI-symmetry in the XY-plane
        IF((PiSymmetry .EQV. .TRUE.) .AND. (DATASETFLAG .GE. 1) .AND. (DATASETFLAG .LE. 3)) THEN
            !Perform 180 degree rotation for vectors
            DO i = minindex(1), (maxindex(1)-1)
                DO j = minindex(2), (maxindex(2)-1)
                    jj = combinedbufsize(1) -i +1
                    kk = combinedbufsize(2) -j +1
                    combinedbuffer(i,j,:) = -1.0D0 * combinedbuffer(jj,kk,:)
                END DO
            END DO
        ELSEIF((PiSymmetry .EQV. .TRUE.) .AND. (DATASETFLAG .LT. 1) .AND. (DATASETFLAG .GT. 3)) THEN
            !Perform 180 degree rotation for non-vectors
            DO i = minindex(1), (maxindex(1)-1)
                DO j = minindex(2), (maxindex(2)-1)
                    jj = combinedbufsize(1) -i +1
                    kk = combinedbufsize(2) -j +1
                    combinedbuffer(i,j,:) = combinedbuffer(jj,kk,:)
                END DO
            END DO
        END IF

        !Z-DIRECTION
        IF( ReflectWhichWay(3) .EQ. +1 ) THEN !DATA IS FOR THE TOP HEMISPHERE
           !The data is now located at the bottom hemisphere, move it to the top hemisphere, remove ghost zones
            DO i = minindex(3), maxindex(3)
                jj = combinedbufsize(3) -i+1
                kk = maxindex(3) + ghost(3) -i+1
                combinedbuffer(:,:,jj) = combinedbuffer(:,:,kk)
            END DO
            !Reflect the top hemisphere values about z = 0
            DO i = minindex(3), (maxindex(3)-1)
                jj = combinedbufsize(3) -i +1
                combinedbuffer(:,:,i) = combinedbuffer(:,:,jj)
            END DO
        ELSE IF( ReflectWhichWay(3) .EQ. -1 ) THEN !DATA IS FOR THE BOTTOM HEMISPHERE
            !Reflect the bottom hemisphere values about z = 0
            DO i = minindex(3), (maxindex(3)-1)
                jj = combinedbufsize(3) -i +1
                combinedbuffer(:,:,jj) = combinedbuffer(:,:,i)
            END DO
        END IF

        !Now, we need to remove the buffer zones and save the results to hdfmetric
        bufzone = ghost*RKsubsteps

        metricdims = combinedbufsize - 2*bufzone

        ALLOCATE( hdfmetric(metricdims(1), metricdims(2), metricdims(3)), STAT=error )
        IF( error .NE. 0 ) STOP "*** Trouble Allocating ***"
        !$OMP PARALLEL DO PRIVATE(j,k,ii,jj,kk)
        DO i = 1, metricdims(1)
            DO j = 1, metricdims(2)
                DO k = 1, metricdims(3)
                    ii = i + bufzone(1)
                    jj = j + bufzone(2)
                    kk = k + bufzone(3)
                    hdfmetric(i,j,k) = combinedbuffer(ii,jj,kk)
                END DO
            END DO
        END DO
        !$OMP END PARALLEL DO

        !Since the extent of grid in the negative direction is equal to the extent of the grid in the positive direction,
        !we can calculate the index for the origin, the value of maximum (xyz), and the value of minimum (xyz) pretty easily.
        Ox = ( metricdims(1) - 1 )/2 + 1
        Oy = ( metricdims(2) - 1 )/2 + 1
        Oz = ( metricdims(3) - 1 )/2 + 1

        Xmax = ( metricdims(1) - Ox ) * delta(1)
        Ymax = ( metricdims(2) - Oy ) * delta(2)
        Zmax = ( metricdims(3) - Oz ) * delta(3)

        Xmin = -1.0D0 * Xmax
        Ymin = -1.0D0 * Ymax
        Zmin = -1.0D0 * Zmax

        RETURN
    END SUBROUTINE ReadHDF5MetricData

