!-----------------------------------------------------------!
! This file consists of six subroutines                     !
! to write and read .dat data for 1,2 and 3-D arrays        !
! (in double precision real number or complex numbers)      !
!-----------------------------------------------------------!

! For reading : the file has to be a long array of 1-D temp array,
! the routine will then transpose the 1-D array into the appropriate dimensions
! For writing: the array is first transposed into 1-D temp array,
! the routine will then write the long 1-D array

! This file also include subroutines to write S and U to file
! with variable filenames that indicate the number of iterations

!========================================================!
!    Read 3D Array Subroutine                            !
!========================================================!

      subroutine Read3d(nx, ny, nz, Array, FileName)
      implicit none

!-----------------------------------------------------------!
!      Declare passed variables                             !
!-----------------------------------------------------------!

         integer*4, intent(in) :: nx, ny, nz
         real*8, intent(out) ::   Array(nx,ny,nz)         
         character*32, intent(in) :: FileName

!-----------------------------------------------------------!
!      Locals                                               !
!-----------------------------------------------------------!

         integer*4   i, j, k
         integer*4   dummy
         integer*4, parameter :: lu = 20 !I/O Unit for Disk I/O
         real*8    temp(nx*ny*nz)

         integer*4 ierror

!-------------------------------------------------------------!
!      Read                                                   !
!-------------------------------------------------------------!

         write(*,1200) FileName
1200     format (' ', 'Reading: ', A, 'now')

         open( UNIT=lu, FILE = FileName, STATUS='OLD', ACTION='READ', &
              &IOSTAT=ierror)

         ! Check to see if OPEN failed
         errorcheck: if (ierror == 6) then
         
            write (*,1020) FileName
1020        format (1X, 'ERROR: File ', A,' does not exist!')
            
         else
                           
            ! File opened successfully, so start reading a 1-D array 'temp'
            readloop: do
               read(lu, *, iostat = ierror) temp        !Get the value
               if ( ierror /= 0 ) exit                  !Exit if not valid
            end do readloop

         end if errorcheck

         dummy = 1

         do i = 1, nx
            do j = 1, ny
               do k = 1, nz
                  Array(i,j,k) = temp(dummy)
                  dummy = dummy + 1
               end do
            end do
         end do

         write (*,1050) FileName
1050     format (' ', 'Reading ', A,' successful')

         return
       end subroutine Read3d


!========================================================!
!    Read 2D Array Subroutine                            !
!========================================================!

       subroutine Read2d(nx, ny, Array, FileName)
        implicit none

!-----------------------------------------------------------!
!      Declare passed variables                             !
!-----------------------------------------------------------!
         
        integer*4, intent(in) :: nx, ny
        real*8, intent(out) ::   Array(nx,ny)
        character*32, intent(in) :: FileName

!-----------------------------------------------------------!
!      Locals                                               !
!-----------------------------------------------------------!

         integer*4   i, j, k
         integer*4   dummy
         integer*4, parameter :: lu = 20 !I/O Unit for Disk I/O
         real*8     temp(nx*ny)

         integer*4 ierror

!-------------------------------------------------------------!
!      Read                                                   !
!-------------------------------------------------------------!

         write(*,1200) FileName
1200     format (' ', 'Reading: ', A, 'now')

         open( UNIT=lu, FILE = FileName, STATUS='OLD', ACTION='READ', &
              &IOSTAT=ierror)

         ! Check to see if OPEN failed
         errorcheck: if (ierror == 6) then
         
            write (*,1020) FileName
1020        format (1X, 'ERROR: File ', A,' does not exist!')
            
         else
                           
         ! File opened successfully, so start reading a 1-D array 'temp'
            readloop: do
               read(lu, *, iostat = ierror) temp        !Get the value
               if ( ierror /= 0 ) exit                  !Exit if not valid
            end do readloop
                                            
         end if errorcheck

         dummy = 1

         do i = 1, nx
            do j = 1, ny
               Array(i,j) = temp(dummy)
               dummy = dummy + 1
            end do
         end do

         write (*,1050) FileName
1050     format (' ', 'Reading ', A,' successful')
         
         return
       end subroutine Read2d


!========================================================!
!    Read 1D Array Subroutine                            !
!========================================================!

      subroutine Read1d(nx, Array, FileName)
        implicit none

!-----------------------------------------------------------!
!      Declare passed variables                             !
!-----------------------------------------------------------!

        integer*4, intent(in) :: nx
        real*8, intent(out) ::   Array(nx)
        character*32, intent(in) :: FileName

!-----------------------------------------------------------!
!      Locals                                               !
!-----------------------------------------------------------!

        integer*4   i, j, k
        integer*4   dummy
        integer*4, parameter :: lu = 20 !I/O Unit for Disk I/O
        real*8     temp(nx)

        integer*4 ierror

!-------------------------------------------------------------!
!      Read                                                   !
!-------------------------------------------------------------!

        write(*,1200) FileName
1200    format (' ', 'Reading: ', A, 'now')

        open( UNIT=lu, FILE = FileName, STATUS='OLD', ACTION='READ', &
             &IOSTAT=ierror)

        ! Check to see if OPEN failed
        errorcheck: if (ierror == 6) then

           write (*,1020) FileName
1020       format (1X, 'ERROR: File ', A,' does not exist!')

        else
                           
           ! File opened successfully, so start reading a 1-D array 'temp'
           readloop: do
              read(lu, *, iostat = ierror) temp !Get the value
              if ( ierror /= 0 ) exit           !Exit if not valid
           end do readloop

        end if errorcheck

        dummy = 1

        do i = 1, nx
           Array(i) = temp(dummy)
           dummy = dummy + 1
        end do

        write (*,1050) FileName
1050    format (' ', 'Reading ', A,' successful')

        return
      end subroutine Read1d

!========================================================!
!    Write 3D Array Subroutine                           !
!========================================================!

      subroutine Write3d(nx, ny, nz, Array, FileName)

        implicit none

!-----------------------------------------------------------!
!      Declare passed variables                             !
!-----------------------------------------------------------!
        
        integer*4, intent(in) :: nx, ny, nz
        real*8, intent(in) ::  Array(nx,ny,nz)
        character*32, intent(in) :: FileName

!-------------------------------------------------------------!
!      Local variables                                        !
!-------------------------------------------------------------!

        integer*4   i, j, k
        integer*4   rank
        integer*4, parameter :: lu = 20 ! I/O unit for disk I/O

        integer*4 ierror        !Iostat variable
        logical :: lexist       !True if file exists
        integer*4 :: lopen      !1 if file is open

!-------------------------------------------------------------!
!      Write                                                  !
!-------------------------------------------------------------!

        write(*,1200) FileName
1200    format (' ', 'Writing: ', A, 'now')
    
        !Does the file already exist?
        lopen = 0

        inquire( FILE=FileName, EXIST=lexist)

        exist:if (.not. lexist) then

           open( UNIT=lu, FILE=FileName, STATUS='NEW', ACTION='WRITE', &
                &IOSTAT=ierror)
           lopen = 1

        else

           !File exists but still replace the file
           open( UNIT=lu, FILE=FileName, STATUS='REPLACE', ACTION='WRITE', &
                &IOSTAT=ierror)
           lopen = 1
        end if exist
		
		
        !Start writing files		
        errorcheck: if (lopen == 0) then

           write(*,1020) FileName
1020       format (' ', 'ERROR: ', A,' cannot be written!')

        else

           do i = 1, nx
              do j = 1, ny
                 do k = 1, nz
                    write(lu, *, iostat=ierror) Array(i,j,k)
                    if( ierror /= 0 ) return            !exit if not valid
                 end do
              end do
           end do

        end if errorcheck

        write(*,1100) FileName
1100    FORMAT(' ', 'Writing ', A, ' successful!')
    
        return
      end subroutine Write3d

!========================================================!
!    Write 2D Array Subroutine                           !
!========================================================!

      subroutine Write2d(nx, ny, Array, FileName)

        implicit none

!-----------------------------------------------------------!
!      Declare passed variables                             !
!-----------------------------------------------------------!

        integer*4, intent(in) :: nx, ny
        real*8, intent(in) ::  Array(nx,ny)
        character*32, intent(in) :: FileName

!-------------------------------------------------------------!
!      Local variables                                        !
!-------------------------------------------------------------!

        integer*4   i, j, k
        integer*4   rank
        integer*4, parameter :: lu = 20 ! I/O unit for disk I/O

        integer*4 ierror        !Iostat variable
        logical :: lexist       !True if file exists
        integer*4 :: lopen      !1 if file is open

!-------------------------------------------------------------!
!      Write                                                  !
!-------------------------------------------------------------!

        write(*,1200) FileName
1200    format (' ', 'Writing: ', A, 'now')

        !Does the file already exist?
        lopen = 0

        inquire( FILE=FileName, EXIST=lexist)

        exist:if (.not. lexist) then

           open( UNIT=lu, FILE=FileName, STATUS='NEW', ACTION='WRITE', &
                &IOSTAT=ierror)
           lopen = 1

        else

           !File exists but still replace the file
           open( UNIT=lu, FILE=FileName, STATUS='REPLACE', ACTION='WRITE', &
                &IOSTAT=ierror)
           lopen = 1
        end if exist


        !Start writing files		
        errorcheck: if (lopen == 0) then

           write(*,1020) FileName
1020       format (' ', 'ERROR: ', A,' cannot be written!')

        else

           do i = 1, nx
              do j = 1, ny
                 write(lu, *, iostat=ierror) Array(i,j)
                 if( ierror /= 0 ) return               !exit if not valid
              end do
           end do


        end if errorcheck

        write(*,1100) FileName
1100    FORMAT(' ', 'Writing ', A, ' successful!')

        return
      end subroutine Write2d
        
!========================================================!
!    Write 1D Array Subroutine                           !
!========================================================!

      subroutine Write1d(nx, Array, FileName)

        implicit none

!-----------------------------------------------------------!
!      Declare passed variables                             !
!-----------------------------------------------------------!

        integer*4, intent(in) :: nx
        real*8, intent(in)      ::  Array(nx)
        character*32, intent(in) :: FileName

!-------------------------------------------------------------!
!      Local variables                                        !
!-------------------------------------------------------------!

        integer*4   i, j, k
        integer*4   rank
        integer*4, parameter :: lu = 20 ! I/O unit for disk I/O

        integer*4 ierror        !Iostat variable
        logical :: lexist       !True if file exists
        integer*4 :: lopen      !1 if file is open

!-------------------------------------------------------------!
!      Write                                                  !
!-------------------------------------------------------------!

        write(*,1200) FileName
1200    format (' ', 'Writing: ', A, 'now')

        !Does the file already exist?
        lopen = 0

        inquire( FILE=FileName, EXIST=lexist)

        exist:if (.not. lexist) then

           open( UNIT=lu, FILE=FileName, STATUS='NEW', ACTION='WRITE', &
                &IOSTAT=ierror)
           lopen = 1

        else

           !File exists but still replace the file
           open( UNIT=lu, FILE=FileName, STATUS='REPLACE', ACTION='WRITE', &
                &IOSTAT=ierror)
           lopen = 1
        end if exist


        !Start writing files		
        errorcheck: if (lopen == 0) then

           write(*,1020) FileName
1020       format (' ', 'ERROR: ', A,' cannot be written!')

        else

           do i = 1, nx
              write(lu, *, iostat=ierror) Array(i)
              if( ierror /= 0 ) return          !exit if not valid
           end do
    
        end if errorcheck

        write(*,1100) FileName
1100    FORMAT(' ', 'Writing ', A, ' successful!')

        return
      end subroutine Write1d

!========================================================!
!    Write 3D Complex Array Subroutine                   !
!========================================================!

      subroutine Write3dC(nx, ny, nz, Array, FileName)

        implicit none

!-----------------------------------------------------------!
!      Declare passed variables                             !
!-----------------------------------------------------------!
        
        integer*4, intent(in) :: nx, ny, nz
        complex*16, intent(in) ::  Array(nx,ny,nz)
        character*32, intent(in) :: FileName

!-------------------------------------------------------------!
!      Local variables                                        !
!-------------------------------------------------------------!

        integer*4   i, j, k
        integer*4   rank
        integer*4, parameter :: lu = 20 ! I/O unit for disk I/O

        integer*4 ierror        !Iostat variable
        logical :: lexist       !True if file exists
        integer*4 :: lopen      !1 if file is open

!-------------------------------------------------------------!
!      Write                                                  !
!-------------------------------------------------------------!

        write(*,1200) FileName
1200    format (' ', 'Writing: ', A, 'now')

        !Does the file already exist?
        lopen = 0

        inquire( FILE=FileName, EXIST=lexist)

        exist:if (.not. lexist) then

           open( UNIT=lu, FILE=FileName, STATUS='NEW', ACTION='WRITE', &
                &IOSTAT=ierror)
           lopen = 1

        else

           !File exists but still replace the file
           open( UNIT=lu, FILE=FileName, STATUS='REPLACE', ACTION='WRITE', &
                &IOSTAT=ierror)
           lopen = 1
        end if exist
		
		
        !Start writing files		
        errorcheck: if (lopen == 0) then

           write(*,1020) FileName
1020       format (' ', 'ERROR: ', A,' cannot be written!')

        else

           do i = 1, nx
              do j = 1, ny
                 do k = 1, nz
                    write(lu, *, iostat=ierror) Array(i,j,k)
                    if( ierror /= 0 ) return            !exit if not valid
                 end do
              end do
           end do
    
    
        end if errorcheck

        write(*,1100) FileName
1100    FORMAT(' ', 'Writing ', A, ' successful!')
    
        return
      end subroutine Write3dC

!========================================================!
!    Write 2D Complex Array Subroutine                   !
!========================================================!

      subroutine Write2dC(nx, ny, Array, FileName)

        implicit none

!-----------------------------------------------------------!
!      Declare passed variables                             !
!-----------------------------------------------------------!

        integer*4, intent(in) :: nx, ny
        complex*16, intent(in) ::  Array(nx,ny)
        character*32, intent(in) :: FileName

!-------------------------------------------------------------!
!      Local variables                                        !
!-------------------------------------------------------------!

        integer*4   i, j, k
        integer*4   rank
        integer*4, parameter :: lu = 20 ! I/O unit for disk I/O

        integer*4 ierror        !Iostat variable
        logical :: lexist       !True if file exists
        integer*4 :: lopen      !1 if file is open

!-------------------------------------------------------------!
!      Write                                                  !
!-------------------------------------------------------------!

        write(*,1200) FileName
1200    format (' ', 'Writing: ', A, 'now')

        !Does the file already exist?
        lopen = 0

        inquire( FILE=FileName, EXIST=lexist)

        exist:if (.not. lexist) then

           open( UNIT=lu, FILE=FileName, STATUS='NEW', ACTION='WRITE', &
                &IOSTAT=ierror)
           lopen = 1

        else
    
           !File exists but still replace the file
           open( UNIT=lu, FILE=FileName, STATUS='REPLACE', ACTION='WRITE', &
                &IOSTAT=ierror)
           lopen = 1
        end if exist


        !Start writing files		
        errorcheck: if (lopen == 0) then

           write(*,1020) FileName
1020       format (' ', 'ERROR: ', A,' cannot be written!')

        else

           do i = 1, nx
              do j = 1, ny
                 write(lu, *, iostat=ierror) Array(i,j)
                 if( ierror /= 0 ) return               !exit if not valid
              end do
           end do


        end if errorcheck

        write(*,1100) FileName
1100    FORMAT(' ', 'Writing ', A, ' successful!')

        return
      end subroutine Write2dC

!========================================================!
!    Write 1D Complex Array Subroutine                   !
!========================================================!

      subroutine Write1dC(nx, Array, FileName)

        implicit none

!-----------------------------------------------------------!
!      Declare passed variables                             !
!-----------------------------------------------------------!

        integer*4, intent(in)    :: nx
        complex*16, intent(in)   :: Array(nx)
        character*32, intent(in) :: FileName

!-------------------------------------------------------------!
!      Local variables                                        !
!-------------------------------------------------------------!

        integer*4   i, j, k
        integer*4   rank
        integer*4, parameter :: lu = 20 ! I/O unit for disk I/O

        integer*4 ierror        !Iostat variable
        logical :: lexist       !True if file exists
        integer*4 :: lopen      !1 if file is open

!-------------------------------------------------------------!
!      Write                                                  !
!-------------------------------------------------------------!

        write(*,1200) FileName
1200    format (' ', 'Writing: ', A, 'now')

        !Does the file already exist?
        lopen = 0

        inquire( FILE=FileName, EXIST=lexist)

        exist:if (.not. lexist) then

           open( UNIT=lu, FILE=FileName, STATUS='NEW', ACTION='WRITE', &
                &IOSTAT=ierror)
           lopen = 1

        else

           !File exists but still replace the file
           open( UNIT=lu, FILE=FileName, STATUS='REPLACE', ACTION='WRITE', &
                &IOSTAT=ierror)
           lopen = 1
        end if exist


        !Start writing files		
        errorcheck: if (lopen == 0) then

           write(*,1020) FileName
1020       format (' ', 'ERROR: ', A,' cannot be written!')

        else

           do i = 1, nx
              write(lu, *, iostat=ierror) Array(i)
              if( ierror /= 0 ) return          !exit if not valid
           end do

        end if errorcheck

        write(*,1100) FileName
1100    FORMAT(' ', 'Writing ', A, ' successful!')

        return
      end subroutine Write1dC

!========================================================!
!    Read 2D Complex Array Subroutine                    !
!========================================================!

       subroutine Read2dC(nx, ny, Array, FileName)
        implicit none

!-----------------------------------------------------------!
!      Declare passed variables                             !
!-----------------------------------------------------------!
         
        integer*4, intent(in) :: nx, ny
        complex*16, intent(out) ::   Array(nx,ny)
        character*32, intent(in) :: FileName

!-----------------------------------------------------------!
!      Locals                                               !
!-----------------------------------------------------------!

         integer*4   i, j, k
         integer*4   dummy
         integer*4, parameter :: lu = 20 !I/O Unit for Disk I/O
         complex*16  temp(nx*ny)

         integer*4 ierror

!-------------------------------------------------------------!
!      Read                                                   !
!-------------------------------------------------------------!

         write(*,1200) FileName
1200     format (' ', 'Reading: ', A, 'now')

         open( UNIT=lu, FILE = FileName, STATUS='OLD', ACTION='READ', &
              &IOSTAT=ierror)
        
         ! Check to see if OPEN failed
         errorcheck: if (ierror == 6) then
         
            write (*,1020) FileName
1020        format (1X, 'ERROR: File ', A,' does not exist!')
        
         else
                           
         ! File opened successfully, so start reading a 1-D array 'temp'
            readloop: do
               read(lu, *, iostat = ierror) temp        !Get the value
               if ( ierror /= 0 ) exit                  !Exit if not valid
            end do readloop
                                            
         end if errorcheck

         dummy = 1

         do i = 1, nx
            do j = 1, ny
               Array(i,j) = temp(dummy)
               dummy = dummy + 1
            end do
         end do

         write (*,1050) FileName
1050     format (' ', 'Reading ', A,' successful')
         
         return
       end subroutine Read2dC

!========================================================!
!    Read 1D Complex Array Subroutine                    !
!========================================================!

      subroutine Read1dC(nx, Array, FileName)
        implicit none

!-----------------------------------------------------------!
!      Declare passed variables                             !
!-----------------------------------------------------------!

        integer*4, intent(in) :: nx
        complex*16, intent(out) ::   Array(nx)
        character*32, intent(in) :: FileName

!-----------------------------------------------------------!
!      Locals                                               !
!-----------------------------------------------------------!

        integer*4   i, j, k
        integer*4   dummy
        integer*4, parameter :: lu = 20 !I/O Unit for Disk I/O
        complex*16  temp(nx)

        integer*4 ierror

!-------------------------------------------------------------!
!      Read                                                   !
!-------------------------------------------------------------!

        write(*,1200) FileName
1200    format (' ', 'Reading: ', A, 'now')

        open( UNIT=lu, FILE = FileName, STATUS='OLD', ACTION='READ', &
             &IOSTAT=ierror)

        ! Check to see if OPEN failed
        errorcheck: if (ierror == 6) then

           write (*,1020) FileName
1020       format (1X, 'ERROR: File ', A,' does not exist!')

        else
                           
           ! File opened successfully, so start reading a 1-D array 'temp'
           readloop: do
              read(lu, *, iostat = ierror) temp !Get the value
              if ( ierror /= 0 ) exit           !Exit if not valid
           end do readloop

        end if errorcheck

        dummy = 1

        do i = 1, nx
           Array(i) = temp(dummy)
           dummy = dummy + 1
        end do

        write (*,1050) FileName
1050    format (' ', 'Reading ', A,' successful')

        return
      end subroutine Read1dC


!========================================================!
!    Write S Array Subroutine                            !
!========================================================!

      SUBROUTINE WriteS(nx, Array,i)

        IMPLICIT none

!-----------------------------------------------------------!
!     Declare passed variables                              !
!-----------------------------------------------------------!

        INTEGER*4, INTENT(IN)::   nx, i
        REAL*8, INTENT(IN)::    Array(nx)

!-----------------------------------------------------------!
!     Declare local variables                               !
!-----------------------------------------------------------!
        
        CHARACTER*32            format_string
        CHARACTER*32            TestFile

!-----------------------------------------------------------!
!     MAIN                                                  !
!-----------------------------------------------------------!

        IF (i < 10) THEN
           format_string = '(A1,I1,A4)'
        ELSE IF (i .ge. 10 .and. i < 100) THEN
           format_string = '(A1,I2,A4)'
        ELSE IF (i < 1000 .and. i .ge. 100) THEN
           format_string = '(A1,I3,A4)'
        ELSE IF (i .ge. 1000 .and. i < 10000) THEN
           format_string = '(A1,I4,A4)'
        ELSE
           format_string = '(A1,I5,A4)'
        END IF

        WRITE(TestFile, format_string) 'S',i,'.dat'

        CALL WRITE1d(nx, Array, TestFile)

        RETURN
      END SUBROUTINE WriteS


!========================================================!
!    Write U Array Subroutine                            !
!========================================================!

      SUBROUTINE WriteU(nx, ny, Array, i)

        IMPLICIT none

!-----------------------------------------------------------!
!     Declare passed variables                              !
!-----------------------------------------------------------!

        INTEGER*4, INTENT(IN)::   nx, ny, i
        REAL*8, INTENT(IN)::    Array(nx,ny)

!-----------------------------------------------------------!
!     Declare local variables                               !
!-----------------------------------------------------------!
        
        CHARACTER*32            format_string
        CHARACTER*32            TestFile

!-----------------------------------------------------------!
!     Write                                                 !
!-----------------------------------------------------------!

        IF (i < 10) THEN
           format_string = '(A1,I1,A4)'
        ELSE IF (i .ge. 10 .and. i < 100) THEN
           format_string = '(A1,I2,A4)'
        ELSE IF (i < 1000 .and. i .ge. 100) THEN
           format_string = '(A1,I3,A4)'
        ELSE IF (i .ge. 1000 .and. i < 10000) THEN
           format_string = '(A1,I4,A4)'
        ELSE
           format_string = '(A1,I5,A4)'
        END IF

        WRITE(TestFile, format_string) 'U',i,'.dat'

        CALL WRITE2d(nx, ny, Array, TestFile)

        RETURN
      END SUBROUTINE WriteU

!========================================================!
!    Write a_nlm Array Subroutine                        !
!========================================================!

      SUBROUTINE Writea(nx, Array, i)

        IMPLICIT none

!-----------------------------------------------------------!
!     Declare passed variables                              !
!-----------------------------------------------------------!

        INTEGER*4, INTENT(IN)::   nx, i
        COMPLEX*16, INTENT(IN)::    Array(nx)

!-----------------------------------------------------------!
!     Declare local variables                               !
!-----------------------------------------------------------!
        
        CHARACTER*32            format_string
        CHARACTER*32            TestFile

!-----------------------------------------------------------!
!     Write                                                 !
!-----------------------------------------------------------!

        IF (i < 10) THEN
           format_string = '(A1,I1,A4)'
        ELSE IF (i .ge. 10 .and. i < 100) THEN
           format_string = '(A1,I2,A4)'
        ELSE IF (i < 1000 .and. i .ge. 100) THEN
           format_string = '(A1,I3,A4)'
        ELSE IF (i .ge. 1000 .and. i < 10000) THEN
           format_string = '(A1,I4,A4)'
        ELSE
           format_string = '(A1,I5,A4)'
        END IF

        WRITE(TestFile, format_string) 'a',i,'.dat'

        CALL WRITE1dC(nx, Array, TestFile)

        RETURN
      END SUBROUTINE Writea

!========================================================!
!    Write Metric data Array Subroutine                  !
!========================================================!

      SUBROUTINE OutputMetric(nx, ny, nz, Array, i, METRICFLAG)

        IMPLICIT none

!-----------------------------------------------------------!
!     Declare passed variables                              !
!-----------------------------------------------------------!

        INTEGER*4, INTENT(IN)::   nx, ny, nz, i, METRICFLAG
        REAL*8, INTENT(IN)::    Array(nx)

!-----------------------------------------------------------!
!     Declare local variables                               !
!-----------------------------------------------------------!
        
        CHARACTER*32            format_string
        CHARACTER*32            TestFile

!-----------------------------------------------------------!
!     Write                                                 !
!-----------------------------------------------------------!

        !METRICFLAG:
        ! = 0 -> alpha, = 1 -> beta1, = 2 -> beta2, = 3 -> beta3
        ! = 4 -> gxx, =5 -> gyy, =6 -> gzz, =7 -> gxy, =8 -> gxz, =9 -> gyz

        IF( METRICFLAG .GE. 0 .AND. METRICFLAG .LE. 3 ) THEN
            IF (i < 10) THEN
                format_string = '(A5,I1,A4)'
            ELSE IF (i .ge. 10 .and. i < 100) THEN
                format_string = '(A5,I2,A4)'
            ELSE IF (i < 1000 .and. i .ge. 100) THEN
                format_string = '(A5,I3,A4)'
            ELSE IF (i .ge. 1000 .and. i < 10000) THEN
                format_string = '(A5,I4,A4)'
            ELSE
                format_string = '(A5,I5,A4)'
            END IF

            IF( METRICFLAG .EQ. 0 ) THEN
                WRITE(TestFile, format_string) 'alpha',i,'.dat'
            ELSEIF( METRICFLAG .EQ. 1 ) THEN
                WRITE(TestFile, format_string) 'beta1',i,'.dat'
            ELSEIF( METRICFLAG .EQ. 2 ) THEN
                WRITE(TestFile, format_string) 'beta2',i,'.dat'
            ELSEIF( METRICFLAG .EQ. 3 ) THEN
                WRITE(TestFile, format_string) 'beta3',i,'.dat'
            END IF

        ELSE
            IF (i < 10) THEN
                format_string = '(A3,I1,A4)'
            ELSE IF (i .ge. 10 .and. i < 100) THEN
                format_string = '(A3,I2,A4)'
            ELSE IF (i < 1000 .and. i .ge. 100) THEN
                format_string = '(A3,I3,A4)'
            ELSE IF (i .ge. 1000 .and. i < 10000) THEN
                format_string = '(A3,I4,A4)'
            ELSE
                format_string = '(A3,I5,A4)'
            END IF

            IF( METRICFLAG .EQ. 4 ) THEN
                WRITE(TestFile, format_string) 'gxx',i,'.dat'
            ELSEIF( METRICFLAG .EQ. 5 ) THEN
                WRITE(TestFile, format_string) 'gyy',i,'.dat'
            ELSEIF( METRICFLAG .EQ. 6 ) THEN
                WRITE(TestFile, format_string) 'gzz',i,'.dat'
            ELSEIF( METRICFLAG .EQ. 7 ) THEN
                WRITE(TestFile, format_string) 'gxy',i,'.dat'
            ELSEIF( METRICFLAG .EQ. 8 ) THEN
                WRITE(TestFile, format_string) 'gxz',i,'.dat'
            ELSEIF( METRICFLAG .EQ. 9 ) THEN
                WRITE(TestFile, format_string) 'gyz',i,'.dat'
            END IF
        END IF

        CALL WRITE3d(nx, ny, nz, Array, TestFile)

        RETURN
      END SUBROUTINE OutputMetric
