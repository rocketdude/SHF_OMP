!----------------------------------------------------------!
!     GetMetric  subroutine                                !
!----------------------------------------------------------!

      SUBROUTINE GetMetric(&
&M, Mr, NP,&
&iter, nchunks,&
&bufsize,&
&r, theta, phi,&
&alpha,&
&betaR, betaTh, betaPhi,&
&gRR, gThTh, gPhiPhi,&
&gRTh, gRPhi, gThPhi)

        USE       omp_lib
        USE       HDF5
        IMPLICIT  none

!-----------------------------------------------------------!
!       Declare calling variables                           !
!-----------------------------------------------------------!

        INTEGER*4               :: M, Mr, NP
        INTEGER*4               :: iter
        INTEGER*4               :: nchunks

        INTEGER(HSIZE_T)        :: bufsize(3)        

        REAL*8                  :: r(Mr+1)
        REAL*8                  :: theta(2*M)
        REAL*8                  :: phi(2*M)

        REAL*8, INTENT(out)     :: alpha(4*NP)
        REAL*8, INTENT(out)     :: betaR(4*NP)
        REAL*8, INTENT(out)     :: betaTh(4*NP)
        REAL*8, INTENT(out)     :: betaPhi(4*NP)

        REAL*8, INTENT(out)     :: gRR(4*NP)
        REAL*8, INTENT(out)     :: gThTh(4*NP)
        REAL*8, INTENT(out)     :: gPhiPhi(4*NP)

        REAL*8, INTENT(out)     :: gRTh(4*NP)
        REAL*8, INTENT(out)     :: gRPhi(4*NP)
        REAL*8, INTENT(out)     :: gThPhi(4*NP)

!-----------------------------------------------------------!
!       Declare local variables                             !
!-----------------------------------------------------------!

        CHARACTER*32            format_string2
        CHARACTER*32            format_string1
        CHARACTER*32            Filename10
        CHARACTER*32            Filename9
        CHARACTER*32            Filename8

        REAL*8                  alphaXYZ(4*NP)
        REAL*8                  betaX(4*NP)
        REAL*8                  betaY(4*NP)
        REAL*8                  betaZ(4*NP)
        REAL*8                  gXX(4*NP)
        REAL*8                  gYY(4*NP)
        REAL*8                  gZZ(4*NP)
        REAL*8                  gXY(4*NP)
        REAL*8                  gXZ(4*NP)
        REAL*8                  gYZ(4*NP)

        REAL*8                  JMatrix(4,4)
        REAL*8                  gcart(4,4)
        REAL*8                  gsph(4,4)
        REAL*8                  dt_ratio

        INTEGER*4               i, j, k
        INTEGER*4               CFLEN2
        INTEGER*4               CFLEN1
        INTEGER*4               crow
        INTEGER*4               error

!----------------------------------------------------------!
!      Main                                                !
!----------------------------------------------------------!

        PRINT *, 'Getting Metric Data of it=', iter

        !The filenames for the metric correspond to the ones defined in the
        !shell program: process-hdf-metric.sh
        
        SELECT CASE(iter)
            CASE( 0:9 )
                CFLEN1 = 9+1+8
                CFLEN2 = 9+1+9
                format_string1 = '(A9,I1,A8)'
                format_string2 = '(A9,I1,A9)'
            CASE( 10:99 )
                CFLEN1 = 9+2+8
                CFLEN2 = 9+2+9
                format_string1 = '(A9,I2,A8)'
                format_string2 = '(A9,I2,A9)'
            CASE( 100:999 )
                CFLEN1 = 9+3+8
                CFLEN2 = 9+3+9
                format_string1 = '(A9,I3,A8)'
                format_string2 = '(A9,I3,A9)'
            CASE( 1000:9999 )
                CFLEN1 = 9+4+8
                CFLEN2 = 9+4+9
                format_string1 = '(A9,I4,A8)'
                format_string2 = '(A9,I4,A9)'
            CASE DEFAULT
                PRINT *, 'Iteration number is too large'
                STOP
        END SELECT

        !ALPHA
        WRITE(Filename9, format_string1) 'alpha.it=',iter,'.rl=9.h5'
        WRITE(Filename8, format_string1) 'alpha.it=',iter,'.rl=8.h5'
        WRITE(Filename10, format_string2) 'alpha.it=',iter,'.rl=10.h5' 
        CALL GetMetricComponent(&
            &M, Mr, NP,&
            &CFLEN2, CFLEN1,&
            &bufsize,&
            &iter, nchunks,&
            &0,&
            &Filename10, Filename9, Filename8,&
            &r, theta, phi,&
            &alphaXYZ)

        !BETA1
        WRITE(Filename9, format_string1) 'beta1.it=',iter,'.rl=9.h5'
        WRITE(Filename8, format_string1) 'beta1.it=',iter,'.rl=8.h5'
        WRITE(Filename10, format_string2) 'beta1.it=',iter,'.rl=10.h5' 
        CALL GetMetricComponent(&
            &M, Mr, NP,&
            &CFLEN2, CFLEN1,&
            &bufsize,&
            &iter, nchunks,&
            &1,&
            &Filename10, Filename9, Filename8,&
            &r, theta, phi,&
            &betaX)

        !BETA2
        WRITE(Filename9, format_string1) 'beta2.it=',iter,'.rl=9.h5'
        WRITE(Filename8, format_string1) 'beta2.it=',iter,'.rl=8.h5'
        WRITE(Filename10, format_string2) 'beta2.it=',iter,'.rl=10.h5' 
        CALL GetMetricComponent(&
            &M, Mr, NP,&
            &CFLEN2, CFLEN1,&
            &bufsize,&
            &iter, nchunks,&
            &2,&
            &Filename10, Filename9, Filename8,&
            &r, theta, phi,&
            &betaY)

        !BETA3
        WRITE(Filename9, format_string1) 'beta3.it=',iter,'.rl=9.h5'
        WRITE(Filename8, format_string1) 'beta3.it=',iter,'.rl=8.h5'
        WRITE(Filename10, format_string2) 'beta3.it=',iter,'.rl=10.h5' 
        CALL GetMetricComponent(&
            &M, Mr, NP,&
            &CFLEN2, CFLEN1,&
            &bufsize,&
            &iter, nchunks,&
            &3,&
            &Filename10, Filename9, Filename8,&
            &r, theta, phi,&
            &betaZ)

        SELECT CASE(iter)
            CASE( 0:9 )
                CFLEN1 = 7+1+8
                CFLEN2 = 7+1+9
                format_string1 = '(A7,I1,A8)'
                format_string2 = '(A7,I1,A9)'
            CASE( 10:99 )
                CFLEN1 = 7+2+8
                CFLEN2 = 7+2+9
                format_string1 = '(A7,I2,A8)'
                format_string2 = '(A7,I2,A9)'
            CASE( 100:999 )
                CFLEN1 = 7+3+8
                CFLEN2 = 7+3+9
                format_string1 = '(A7,I3,A8)'
                format_string2 = '(A7,I3,A9)'
            CASE( 1000:9999 )
                CFLEN1 = 7+4+8
                CFLEN2 = 7+4+9
                format_string1 = '(A7,I4,A8)'
                format_string2 = '(A7,I4,A9)'
            CASE DEFAULT
                PRINT *, 'Iteration number is too large'
                STOP
        END SELECT

        !GXX
        WRITE(Filename9, format_string1) 'gxx.it=',iter,'.rl=9.h5'
        WRITE(Filename8, format_string1) 'gxx.it=',iter,'.rl=8.h5'
        WRITE(Filename10, format_string2) 'gxx.it=',iter,'.rl=10.h5' 
        CALL GetMetricComponent(&
            &M, Mr, NP,&
            &CFLEN2, CFLEN1,&
            &bufsize,&
            &iter, nchunks,&
            &4,&
            &Filename10, Filename9, Filename8,&
            &r, theta, phi,&
            &gXX)

        !GYY
        WRITE(Filename9, format_string1) 'gyy.it=',iter,'.rl=9.h5'
        WRITE(Filename8, format_string1) 'gyy.it=',iter,'.rl=8.h5'
        WRITE(Filename10, format_string2) 'gyy.it=',iter,'.rl=10.h5' 
        CALL GetMetricComponent(&
            &M, Mr, NP,&
            &CFLEN2, CFLEN1,&
            &bufsize,&
            &iter, nchunks,&
            &5,&
            &Filename10, Filename9, Filename8,&
            &r, theta, phi,&
            &gYY)

        !GZZ
        WRITE(Filename9, format_string1) 'gzz.it=',iter,'.rl=9.h5'
        WRITE(Filename8, format_string1) 'gzz.it=',iter,'.rl=8.h5'
        WRITE(Filename10, format_string2) 'gzz.it=',iter,'.rl=10.h5' 
        CALL GetMetricComponent(&
            &M, Mr, NP,&
            &CFLEN2, CFLEN1,&
            &bufsize,&
            &iter, nchunks,&
            &6,&
            &Filename10, Filename9, Filename8,&
            &r, theta, phi,&
            &gZZ)

        !GXY
        WRITE(Filename9, format_string1) 'gxy.it=',iter,'.rl=9.h5'
        WRITE(Filename8, format_string1) 'gxy.it=',iter,'.rl=8.h5'
        WRITE(Filename10, format_string2) 'gxy.it=',iter,'.rl=10.h5' 
        CALL GetMetricComponent(&
            &M, Mr, NP,&
            &CFLEN2, CFLEN1,&
            &bufsize,&
            &iter, nchunks,&
            &7,&
            &Filename10, Filename9, Filename8,&
            &r, theta, phi,&
            &gXY)

        !GXZ
        WRITE(Filename9, format_string1) 'gxz.it=',iter,'.rl=9.h5'
        WRITE(Filename8, format_string1) 'gxz.it=',iter,'.rl=8.h5'
        WRITE(Filename10, format_string2) 'gxz.it=',iter,'.rl=10.h5' 
        CALL GetMetricComponent(&
            &M, Mr, NP,&
            &CFLEN2, CFLEN1,&
            &bufsize,&
            &iter, nchunks,&
            &8,&
            &Filename10, Filename9, Filename8,&
            &r, theta, phi,&
            &gXZ)

        !GYZ
        WRITE(Filename9, format_string1) 'gyz.it=',iter,'.rl=9.h5'
        WRITE(Filename8, format_string1) 'gyz.it=',iter,'.rl=8.h5'
        WRITE(Filename10, format_string2) 'gyz.it=',iter,'.rl=10.h5' 
        CALL GetMetricComponent(&
            &M, Mr, NP,&
            &CFLEN2, CFLEN1,&
            &bufsize,&
            &iter, nchunks,&
            &9,&
            &Filename10, Filename9, Filename8,&
            &r, theta, phi,&
            &gYZ)

        !Now we need to perform coordinate transformation from (t,x,y,z) to (t,r,th,phi)
        !Both g_sph and g_cart are of down indices
 
        !$OMP PARALLEL DO PRIVATE(j, k, JMatrix, crow, gcart, gsph, error)
        DO i = 1, (Mr+1)        
            DO j = 1, (2*M)
                DO k = 1, (2*M)
            
                    !Build the Jacobian matrix JMatrix
                    CALL EvaluateJacobian(r(i),theta(j),phi(k),JMatrix)

                    crow = (i-1)*(2*M)*(2*M) + (j-1)*(2*M) + k
                   
                    ! Build the matrix g_cart
                    CALL EvaluateMatrixofMetric(&
                         &alphaXYZ(crow), betaX(crow), betaY(crow), betaZ(crow),&
                         &gXX(crow), gYY(crow), gZZ(crow),&
                         &gXY(crow), gXZ(crow), gYZ(crow),&
                         &gcart)

                    !  g_sph = (J^T) * g_cart * J
                    gsph = MATMUL(TRANSPOSE(JMatrix), MATMUL(gcart,JMatrix))

                    ! Get the 3-metric by inverting the 3-matric part of g_sph
                    CALL Invert3Metric(gsph(2,2), gsph(3,3), gsph(4,4),&
                                      &gsph(2,3), gsph(2,4), gsph(3,4),&
                                      &gRR(crow), gThTh(crow), gPhiPhi(crow),&
                                      &gRTh(crow), gRPhi(crow), gThPhi(crow),&
                                      &error)

                    IF( error .NE. 0 ) THEN
                        PRINT *, 'Determinant is zero, error in inverting metric'
                        PRINT *, 'i=',i,',j=',j,',k=',k
                        STOP
                    END IF

                    ! Get beta^r, beta^th, beta^phi
                    ! gsph(1,2) = g_01 = beta_r, gsph(1,3) = g_02 = beta_th, gsph(1,4) = g_03 = beta_phi
                    betaR(crow) = gRR(crow)*gsph(1,2) + gRTh(crow)*gsph(1,3) + gRPhi(crow)*gsph(1,4)
                    betaTh(crow) = gRTh(crow)*gsph(1,2) + gThTh(crow)*gsph(1,3) + gThPhi(crow)*gsph(1,4)
                    betaPhi(crow) = gRPhi(crow)*gsph(1,2) + gThPhi(crow)*gsph(1,3) + gPhiPhi(crow)*gsph(1,4)

                    ! Get alpha
                    ! alpha = sqrt( beta^i beta_i - g_00 )
                    ! beta^r = betaR, beta^th = betaTh, beta^phi = betaPhi, gsph(1,1) = g_00
                    alpha(crow) = SQRT( gsph(1,2) * betaR(crow) + gsph(1,3) * betaTh(crow) + gsph(1,4) * betaPhi(crow) &
                                & - gsph(1,1) )

                END DO
            END DO
        END DO
        !$OMP END PARALLEL DO

        PRINT *, 'DONE!'

        RETURN
      END SUBROUTINE GetMetric

