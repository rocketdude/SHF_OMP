!----------------------------------------------------------!
!     GetMetric  subroutine                                !
!----------------------------------------------------------!

      SUBROUTINE GetMetric(&
&M, Mr, NP,&
&h10, h9,&
&iter, nchunks,&
&bufsize,&
&r, theta, phi,&
&alpha,&
&betaR, betaTh, betaPhi,&
&gRR, gThTh, gPhiPhi,&
&gRTh, gRPhi, gThPhi)

        USE       omp_lib
        IMPLICIT  none

!-----------------------------------------------------------!
!       Declare calling variables                           !
!-----------------------------------------------------------!

        INTEGER*4, INTENT(in)  :: M, Mr, NP
        INTEGER*4, INTENT(in)  :: iter
        INTEGER*4, INTENT(in)  :: nchunks

        INTEGER(HSIZE_T), INTENT(in):: bufsize(3)        

        REAL*8, INTENT(in)   :: h10, h9
        REAL*8, INTENT(in)   :: r(Mr+1)
        REAL*8, INTENT(in)   :: theta(2*M)
        REAL*8, INTENT(in)   :: phi(2*M)

        REAL*8, INTENT(out)  :: alpha(4*NP)
        REAL*8, INTENT(out)  :: betaR(4*NP)
        REAL*8, INTENT(out)  :: betaTh(4*NP)
        REAL*8, INTENT(out)  :: betaPhi(4*NP)

        REAL*8, INTENT(out)  :: gRR(4*NP)
        REAL*8, INTENT(out)  :: gThTh(4*NP)
        REAL*8, INTENT(out)  :: gPhiPhi(4*NP)

        REAL*8, INTENT(out)  :: gRTh(4*NP)
        REAL*8, INTENT(out)  :: gRPhi(4*NP)
        REAL*8, INTENT(out)  :: gThPhi(4*NP)

!-----------------------------------------------------------!
!       Declare local variables                             !
!-----------------------------------------------------------!

        INTEGER*4               LM(64, 64)

        CHARACTER*32            format_string10
        CHARACTER*32            format_string9
        CHARACTER*32            Filename10
        CHARACTER*32            Filename9

        REAL*8                  Calpha(4*NP)
        REAL*8                  betaX(4*NP)
        REAL*8                  betaY(4*NP)
        REAL*8                  betaZ(4*NP)
        REAL*8                  gXX(4*NP)
        REAL*8                  gYY(4*NP)
        REAL*8                  gZZ(4*NP)
        REAL*8                  gXY(4*NP)
        REAL*8                  gXZ(4*NP)
        REAL*8                  gYZ(4*NP)

        REAL*8                  TMatrix(4,4)
        REAL*8                  gcart(4,4)
        REAL*8                  gsph(4,4)

!----------------------------------------------------------!
!      Main                                                !
!----------------------------------------------------------!

        PRINT *, 'Getting Metric Data'

        CALL GetLekienCoefficients(LM)

        !The filenames for the metric correspond to the ones defined in the
        !shell program: process-hdf-metric.sh
        
        SELECT CASE(iter)
            CASE( 0:9 )
                format_string9 = '(A9,I1,A7)'
                format_string10 = '(A9,I1,A8)'
            CASE( 10:99 )
                format_string9 = '(A9,I2,A7)'
                format_string10 = '(A9,I2,A8)'
            CASE( 100:999 )
                format_string9 = '(A9,I3,A7)'
                format_string10 = '(A9,I3,A8)'
            CASE( 1000:9999 )
                format_string9 = '(A9,I3,A7)'
                format_string10 = '(A9,I3,A8)'
            CASE DEFAULT
                PRINT *, 'Iteration number is too large'
                STOP
        END SELECT

        !ALPHA
        WRITE(Filename9, format_string9) 'alpha.it=',iter,'rl=9.h5'
        WRITE(Filename10, format_string10) 'alpha.it=',iter,'rl=10.h5' 
        CALL ReadNInterpolateMetric(&
            &M, Mr, NP,&
            &LM,&
            &bufsize,&
            &iter, nchunks,&
            &0,&
            &Filename10, Filename9,&
            &r, theta, phi,&
            &Calpha)

        !BETA1
        WRITE(Filename9, format_string9) 'beta1.it=',iter,'rl=9.h5'
        WRITE(Filename10, format_string10) 'beta1.it=',iter,'rl=10.h5' 
        CALL ReadNInterpolateMetric(&
            &M, Mr, NP,&
            &LM,&
            &bufsize,&
            &iter, nchunks,&
            &1,&
            &Filename10, Filename9,&
            &r, theta, phi,&
            &betaX)

        !BETA2
        WRITE(Filename9, format_string9) 'beta2.it=',iter,'rl=9.h5'
        WRITE(Filename10, format_string10) 'beta2.it=',iter,'rl=10.h5' 
        CALL ReadNInterpolateMetric(&
            &M, Mr, NP,&
            &LM,&
            &bufsize,&
            &iter, nchunks,&
            &2,&
            &Filename10, Filename9,&
            &r, theta, phi,&
            &betaY)

        !BETA3
        WRITE(Filename9, format_string9) 'beta3.it=',iter,'rl=9.h5'
        WRITE(Filename10, format_string10) 'beta3.it=',iter,'rl=10.h5' 
        CALL ReadNInterpolateMetric(&
            &M, Mr, NP,&
            &LM,&
            &bufsize,&
            &iter, nchunks,&
            &3,&
            &Filename10, Filename9,&
            &r, theta, phi,&
            &betaZ)

        SELECT CASE(iter)
            CASE( 0:9 )
                format_string9 = '(A7,I1,A7)'
                format_string10 = '(A7,I1,A8)'
            CASE( 10:99 )
                format_string9 = '(A7,I2,A7)'
                format_string10 = '(A7,I2,A8)'
            CASE( 100:999 )
                format_string9 = '(A7,I3,A7)'
                format_string10 = '(A7,I3,A8)'
            CASE( 1000:9999 )
                format_string9 = '(A7,I3,A7)'
                format_string10 = '(A7,I3,A8)'
            CASE DEFAULT
                PRINT *, 'Iteration number is too large'
                STOP
        END SELECT

        !GXX
        WRITE(Filename9, format_string9) 'gxx.it=',iter,'rl=9.h5'
        WRITE(Filename10, format_string10) 'gxx.it=',iter,'rl=10.h5' 
        CALL ReadNInterpolateMetric(&
            &M, Mr, NP,&
            &LM,&
            &bufsize,&
            &iter, nchunks,&
            &4,&
            &Filename10, Filename9,&
            &r, theta, phi,&
            &gXX)

        !GYY
        WRITE(Filename9, format_string9) 'gyy.it=',iter,'rl=9.h5'
        WRITE(Filename10, format_string10) 'gyy.it=',iter,'rl=10.h5' 
        CALL ReadNInterpolateMetric(&
            &M, Mr, NP,&
            &LM,&
            &bufsize,&
            &iter, nchunks,&
            &5,&
            &Filename10, Filename9,&
            &r, theta, phi,&
            &gYY)

        !GZZ
        WRITE(Filename9, format_string9) 'gzz.it=',iter,'rl=9.h5'
        WRITE(Filename10, format_string10) 'gzz.it=',iter,'rl=10.h5' 
        CALL ReadNInterpolateMetric(&
            &M, Mr, NP,&
            &LM,&
            &bufsize,&
            &iter, nchunks,&
            &6,&
            &Filename10, Filename9,&
            &r, theta, phi,&
            &gZZ)

        !GXY
        WRITE(Filename9, format_string9) 'gxy.it=',iter,'rl=9.h5'
        WRITE(Filename10, format_string10) 'gxy.it=',iter,'rl=10.h5' 
        CALL ReadNInterpolateMetric(&
            &M, Mr, NP,&
            &LM,&
            &bufsize,&
            &iter, nchunks,&
            &7,&
            &Filename10, Filename9,&
            &r, theta, phi,&
            &gXY)

        !GXZ
        WRITE(Filename9, format_string9) 'gxz.it=',iter,'rl=9.h5'
        WRITE(Filename10, format_string10) 'gxz.it=',iter,'rl=10.h5' 
        CALL ReadNInterpolateMetric(&
            &M, Mr, NP,&
            &LM,&
            &bufsize,&
            &iter, nchunks,&
            &8,&
            &Filename10, Filename9,&
            &r, theta, phi,&
            &gXZ)

        !GYZ
        WRITE(Filename9, format_string9) 'gyz.it=',iter,'rl=9.h5'
        WRITE(Filename10, format_string10) 'gyz.it=',iter,'rl=10.h5' 
        CALL ReadNInterpolateMetric(&
            &M, Mr, NP,&
            &LM,&
            &bufsize,&
            &iter, nchunks,&
            &9,&
            &Filename10, Filename9,&
            &r, theta, phi,&
            &gYZ)

        !Now we need to perform coordinate transformation from (t,x,y,z) to (t,r,th,phi)

        DO j = 1, (2*M)
            DO k = 1, (2*M)

                !Build the transformation matrix: TMatrix (independent of coordinate r)
                CALL EvaluateTransformationMatrix(theta(j),phi(k),TMatrix)

                DO i = 1, (Mr+1)
            
                    crow = (i-1)*(2*M)*(2*M) + (j-1)*(2*M) + k

                    
                    

        PRINT *, 'DONE!'

        RETURN
      END SUBROUTINE GetMetric
