!----------------------------------------------------------!
!     ReadHDF5Metric  subroutine                           !
!----------------------------------------------------------!

    SUBROUTINE ReadHDF5Metric(&
        &Filename, dataset,&
        &nchunks,&
        &bufsize,&
        &Xmin, Ymin, Zmin,&
        &Xmax, Ymax, Zmax,&
        &delta,&
        &Ox, Oy, Oz,&
        &metricdims,&
        &hdfmetric)


        !This subroutine is to read the metric data supplied in HDF5 file.
        !Since the metric data is chunked into pieces, the subroutine
        !also recombines all the data into a single 3-dimensional array hdfmetric.
        !Our algorithm of recombining the metric data removes the outer ghost zones
        !whenever we have to reflect the data.

        USE omp_lib
        USE HDF5
        USE h5lt
        IMPLICIT  none

!-----------------------------------------------------------!
!       Declare calling variables                           !
!-----------------------------------------------------------!

        INTEGER*4, INTENT(in)               :: nchunks
        INTEGER(HSIZE_T), INTENT(in)        :: bufsize(3)

        CHARACTER*32, INTENT(in)            :: Filename
        CHARACTER*32, INTENT(in)            :: dataset(nchunks)

        REAL*8, INTENT(out)                 :: Xmin, Ymin, Zmin
        REAL*8, INTENT(out)                 :: Xmax, Ymax, Zmax
        REAL*8, INTENT(out)                 :: delta(3)
        INTEGER*4, INTENT(out)              :: metricdims(3)
        REAL*8, DIMENSION(:,:,:), POINTER   :: hdfmetric

!-----------------------------------------------------------!
!       Declare local variables                             !
!-----------------------------------------------------------!

        !Used to read the HDF5 files
        INTEGER(HID_T) ::           file_id
        INTEGER(HID_T) ::           dims(3) !The dimensions of a dataset
        INTEGER(HID_T) ::           nxyz(nchunks,3) !The dimensions of the dataset for each chunk
        INTEGER*4                   hdferr  !Lets us know if there is error in reading the hdf5 files

        REAL*8                      buffer(nchunks,bufsize(1),bufsize(2),bufsize(3))
                                    !Buffer is the buffer for reading the dataset: buffer(chunk #, x, y, z)
                                    !Increase the buffer size if the dataset is larger, otherwise data will be truncated

        !Used to chunk together the chunked pieces
        INTEGER*4                   ReflectWhichWay(3)  ! 0 --> no need to reflect
                                                        ! -1 --> reflect the bottom (negative) data to the top
                                                        ! +1 --> reflect the top (positive) data to the bottom
        INTEGER*4                   iorigin(nchunks,3)  !Bottom left index for each dataset
        INTEGER*4                   minindex(3)
        INTEGER*4                   maxindex(3)
        INTEGER*4                   ghost(3)            !Extra grids due to ghost zone
        REAL*8                      origin(nchunks,3)   !Bottom left values for each dataset
        REAL*8                      minpoint(3)
        REAL*8                      maxpoint(3)
        REAL*8                      fudge               !fudge factor

        !just placeholders in reading the HDF5 files--not really used
        INTEGER*4       ::          type_class
        INTEGER(SIZE_T) ::          type_size

        !DO loop counters
        INTEGER*4                   i, j, k
        INTEGER*4                   ii, jj, kk
        INTEGER*4                   cnum                !chunk number

!----------------------------------------------------------!
!      Main                                                !
!----------------------------------------------------------!

        fudge = 1.0D-7

        !Read the HDF5 files
        hdferr = 0

        !Open the file
        CALL h5fopen_f(Filename, H5F_ACC_RDONLY_F, file_id, hdferr, H5P_DEFAULT_F)
        IF( hdferr .NE. 0 ) THEN
            PRINT *, 'ERROR in opening HDF5 file'
            STOP
        END IF

        !Iterate through all the chunks
        DO cnum = 1, nchunks

            CALL h5ltread_dataset_double_f(file_id, dataset(cnum), H5T_IEEE_F64LE, bufsize, buffer(cnum,:,:,:), hdferr) 
            IF( hdferr .NE. 0 ) THEN
                PRINT *, 'ERROR in reading dataset'
                STOP
            END IF
            CALL h5ltget_dataset_info_f(file_id, dataset(cnum), dims, type_class, type_size, hdferr)          
            IF( hdferr .NE. 0 ) THEN
                PRINT *, 'ERROR in getting the dimensions of the dataset'
                STOP
            END IF
            CALL h5ltget_attribute_int_f(file_id, dataset(cnum), 'iorigin', hdferr, iorigin(cnum,:)) 
            IF( hdferr .NE. 0 ) THEN
                PRINT *, 'ERROR in getting the iorigin'
                STOP
            END IF
            CALL h5ltget_attribute_double_f(file_id, dataset(cnum), 'origin', hdferr, origin(cnum,:)) 
            IF( hdferr .NE. 0 ) THEN
                PRINT *, 'ERROR in getting the values of the origin'
                STOP
            END IF

            nxyz(cnum,:) = dims(:)
            IF( (dims(1) .GT. bufsize(1)) .OR. (dims(2) .GT. bufsize(2)) .OR. (dims(3) .GT. bufsize(3)) ) THEN
                PRINT *, 'ERROR: buffer size is not big enough'
                STOP
            END IF

        END DO

        !All the chunks should have the same spatial discretizations dx, dy, dz, so just get these delta from a single dataset
        CALL h5ltget_attribute_double_f(file_id, dataset(nchunks), 'delta', hdferr, delta)
        IF( hdferr .NE. 0 ) THEN
            PRINT *, 'ERROR in getting dx, dy, dz'
            STOP
        END IF
        !Close the file
        CALL h5fclose_f(file_id, hdferr)
        IF( hdferr .NE. 0 ) THEN
            PRINT *, 'ERROR in closing the HDF5 file'
            STOP
        END IF

        DO i = 1,3
            !Find the left-bottom corner (minimum index)
            minindex(i) = MINVAL(iorigin(:,i))

            !Find the top-right corner (maximum index)
            maxindex(i) = MAXVAL(iorigin(:,i)+nxyz(:,i))

            !Make minindex = (1,1,1) (left-bottom) and change the maxindex & iorigin accordingly
            iorigin(:,i) = iorigin(:,i) - minindex(i) +1
            maxindex(i) = maxindex(i) - minindex(i) +1
            minindex(i) = 1

            !Find the minimum x,y,and z
            minpoint(i) = MINVAL(origin(:,i))
            !Find the maximum x,y,and z
            maxpoint(i) = minpoint(i) + DBLE(maxindex(i)-minindex(i))*delta(i)

            !Find out if we either only have half the data or all the data, in each direction
            IF( ABS(maxpoint(i)) .EQ. ABS(minpoint(i)) ) THEN
                metricdims(i) = maxindex(i)
                ReflectWhichWay(i) = 0
            ELSE IF( ABS(maxpoint(i)) .GT. ABS(minpoint(i)) ) THEN
                IF( minpoint(i) .GT. 0 ) THEN
                    PRINT *, 'ERROR: in spatial index', i, ', the minimum point of the metric has to be less than 0'
                    STOP
                ELSE
                    ghost(i) = IDINT(fudge + ABS(minpoint(i))/delta(i) ) !Any data below the 0 plane would be ghost zone
                    maxindex(i) = maxindex(i) - ghost(i)    !remove the bottom ghost zones but keep minindex = 1
                    metricdims(i) = 2*( maxindex(i)-minindex(i) ) + 1
                    ReflectWhichWay(i) = +1
                END IF
            ELSE IF( ABS(maxpoint(i)) .LT. ABS(minpoint(i)) ) THEN
                IF( maxpoint(i) .LT. 0 ) THEN
                    PRINT *, 'ERROR: in spatial index', i, ', the maximum point of the metric has to be more than 0'
                    STOP
                ELSE
                    ghost(i) = IDINT(fudge + maxpoint(i)/delta(i))  !Any data above the 0 plane would be ghost zone
                    maxindex(i) = maxindex(i) - ghost(i)     !remove the top ghost zones and keep minindex = 1
                    metricdims(i) = 2*( maxindex(i)-minindex(i) ) + 1 
                    ReflectWhichWay(i) = -1
                END IF
            END IF
        END DO

        !Allocate hdfmetric to its appropriate size
        ALLOCATE( hdfmetric(metricdims(1), metricdims(2), metricdims(3)) )

        !Combine the chunked data and save it in hdfmetric
        DO cnum = 1, nchunks
            DO i = 1, nxyz(cnum,1)
                DO j = 1, nxyz(cnum,2)
                    DO k = 1, nxyz(cnum,3)
                        
                        hdfmetric( (i+iorigin(cnum,1)-1),&
                                &(j+iorigin(cnum,2)-1),&
                                &(k+iorigin(cnum,3)-1) ) = buffer(cnum,i,j,k)   !NOTE:hdfmetric still includes ghost zones

                    END DO
                END DO
            END DO
        END DO

        !Perform the reflection accordingly
        !NOTE: There is no direct relationship between (i,j,k)/(ii,jj,kk) and the coordinates (x,y,z)
        !(i,j,k) and (ii,jj,kk) are really used as dummy place holders here
        
        !X-DIRECTION
        IF( ReflectWhichWay(1) .EQ. +1 ) THEN !DATA IS FOR THE TOP HEMISPHERE
            !Transfer the values of the metric (minus the ghost zones) to the top
            !metricdims - maxindex +1 should correspond to x = 0
            DO i = minindex(1), maxindex(1)
                ii = maxindex(1) +  ghost(1) -i +1    !Start from the higest index as if the ghost zones are still there
                jj = metricdims(1) -i +1
                hdfmetric(jj,:,:) = hdfmetric(ii,:,:)
            END DO
            !The values of the metric below x = 0 is worthless now, reflect the values about x = 0
            DO i = minindex(1), (maxindex(1)-1)
                jj = metricdims(1) -i +1
                hdfmetric(i,:,:) = hdfmetric(jj,:,:)
            END DO
        ELSE IF( ReflectWhichWay(1) .EQ. -1 ) THEN !DATA IS FOR THE BOTTOM HEMISPHERE
            !The values of the metric above x = 0 is worthless, reflect the values about x = 0
            DO i = minindex(1), (maxindex(1)-1)
                jj = metricdims(1) -i +1
                hdfmetric(jj,:,:) = hdfmetric(i,:,:)
            END DO
        END IF

        !Y-DIRECTION
        IF( ReflectWhichWay(2) .EQ. +1 ) THEN !DATA IS FOR THE TOP HEMISPHERE
            !Transfer the values of the metric (minus the ghost zones) to the top
            !metricdims - maxindex +1 should correspond to y = 0
            DO i = minindex(2), maxindex(2)
                ii = maxindex(2) +  ghost(2) -i +1    !Start from the higest index as if the ghost zones are still there
                jj = metricdims(2) -i +1
                hdfmetric(:,jj,:) = hdfmetric(:,ii,:)
            END DO
            !The values of the metric below y = 0 is worthless now, reflect the values about y = 0
            DO i = minindex(2), (maxindex(2)-1)
                jj = metricdims(2) -i +1
                hdfmetric(:,i,:) = hdfmetric(:,jj,:)
            END DO
        ELSE IF( ReflectWhichWay(2) .EQ. -1 ) THEN !DATA IS FOR THE BOTTOM HEMISPHERE
            !The values of the metric above y = 0 is worthless, reflect the values about y = 0
            DO i = minindex(2), (maxindex(2)-1)
                jj = metricdims(2) -i +1
                hdfmetric(:,jj,:) = hdfmetric(:,i,:)
            END DO
        END IF
        
        !Z-DIRECTION
        IF( ReflectWhichWay(3) .EQ. +1 ) THEN !DATA IS FOR THE TOP HEMISPHERE
            !Transfer the values of the metric (minus the ghost zones) to the top
            !metricdims - maxindex +1 should correspond to z = 0
            DO i = minindex(3), maxindex(3)
                ii = maxindex(3) +  ghost(3) -i +1    !Start from the higest index as if the ghost zones are still there
                jj = metricdims(3) -i +1
                hdfmetric(:,:,jj) = hdfmetric(:,:,ii)
            END DO
            !The values of the metric below z = 0 is worthless now, reflect the values about z = 0
            DO i = minindex(3), (maxindex(3)-1)
                jj = metricdims(3) -i +1
                hdfmetric(:,:,i) = hdfmetric(:,:,jj)
            END DO
        ELSE IF( ReflectWhichWay(3) .EQ. -1 ) THEN !DATA IS FOR THE BOTTOM HEMISPHERE
            !The values of the metric above z = 0 is worthless, reflect the values about z = 0
            DO i = minindex(3), (maxindex(3)-1)
                jj = metricdims(3) -i +1
                hdfmetric(:,:,jj) = hdfmetric(:,:,jj)
            END DO
        END IF

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

    END SUBROUTINE ReadHDF5Metric
