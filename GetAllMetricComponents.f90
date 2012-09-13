!----------------------------------------------------------!
!     GetMetric  subroutine                                !
!----------------------------------------------------------!

      SUBROUTINE GetAllMetricComponents(&
&M, Mr, NP,&
&iter, nchunks,&
&bufsize,&
&time,&
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

        REAL*8, INTENT(out)     :: time

!-----------------------------------------------------------!
!       Declare local variables                             !
!-----------------------------------------------------------!

        CHARACTER*32            format_string_LapseShift(nchunks)
        CHARACTER*32            format_string_SpatMetric(nchunks)

        CHARACTER*32            alphafilename(nchunks)
        CHARACTER*32            beta1filename(nchunks)
        CHARACTER*32            beta2filename(nchunks)
        CHARACTER*32            beta3filename(nchunks)
        CHARACTER*32            gxxfilename(nchunks)
        CHARACTER*32            gyyfilename(nchunks)
        CHARACTER*32            gzzfilename(nchunks)
        CHARACTER*32            gxyfilename(nchunks)
        CHARACTER*32            gxzfilename(nchunks)
        CHARACTER*32            gyzfilename(nchunks)

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
        REAL*8                  timeofdata(10)

        INTEGER*4               i, j, k
        INTEGER*4               crow
        INTEGER*4               error

!----------------------------------------------------------!
!      Main                                                !
!----------------------------------------------------------!

        PRINT *, 'Getting Metric Data of it=', iter

        !$OMP PARALLEL DO PRIVATE(j)
        DO i = 1,nchunks
            j = i - 1
            IF( j .LT. 10 ) THEN
                format_string_LapseShift(i) = '(A11,I1,A3)'
                format_string_SpatMetric(i) = '(A9,I1,A3)'
            ELSE
                format_string_LapseShift(i) = '(A11,I2,A3)'
                format_string_SpatMetric(i) = '(A9,I2,A3)'
            END IF
            
            WRITE(alphafilename(i), format_string_LapseShift(i)) 'alpha.file_', j, '.h5'
            WRITE(beta1filename(i), format_string_LapseShift(i)) 'beta1.file_', j, '.h5'
            WRITE(beta2filename(i), format_string_LapseShift(i)) 'beta2.file_', j, '.h5'
            WRITE(beta3filename(i), format_string_LapseShift(i)) 'beta3.file_', j, '.h5'
            WRITE(gxxfilename(i), format_string_SpatMetric(i)) 'gxx.file_', j, '.h5'
            WRITE(gyyfilename(i), format_string_SpatMetric(i)) 'gyy.file_', j, '.h5'
            WRITE(gzzfilename(i), format_string_SpatMetric(i)) 'gzz.file_', j, '.h5'
            WRITE(gxyfilename(i), format_string_SpatMetric(i)) 'gxy.file_', j, '.h5'
            WRITE(gxzfilename(i), format_string_SpatMetric(i)) 'gxz.file_', j, '.h5'
            WRITE(gyzfilename(i), format_string_SpatMetric(i)) 'gyz.file_', j, '.h5'
        END DO
        !$OMP END PARALLEL DO
                
        !ALPHA
        CALL GetMetricComponent(&
            &M, Mr, NP,&
            &bufsize,&
            &timeofdata(1),&
            &iter, nchunks,&
            &0,&
            &alphafilename,&
            &r, theta, phi,&
            &alphaXYZ)

        !BETA1
        CALL GetMetricComponent(&
            &M, Mr, NP,&
            &bufsize,&
            &timeofdata(2),&
            &iter, nchunks,&
            &1,&
            &beta1filename,&
            &r, theta, phi,&
            &betaX)

        !BETA2
        CALL GetMetricComponent(&
            &M, Mr, NP,&
            &bufsize,&
            &timeofdata(3),&
            &iter, nchunks,&
            &2,&
            &beta2filename,&
            &r, theta, phi,&
            &betaY)

        !BETA3
        CALL GetMetricComponent(&
            &M, Mr, NP,&
            &bufsize,&
            &timeofdata(4),&
            &iter, nchunks,&
            &3,&
            &beta3filename,&
            &r, theta, phi,&
            &betaZ)

        !GXX
        CALL GetMetricComponent(&
            &M, Mr, NP,&
            &bufsize,&
            &timeofdata(5),&
            &iter, nchunks,&
            &4,&
            &gxxfilename,&
            &r, theta, phi,&
            &gXX)

        !GYY
        CALL GetMetricComponent(&
            &M, Mr, NP,&
            &bufsize,&
            &timeofdata(6),&
            &iter, nchunks,&
            &5,&
            &gyyfilename,&
            &r, theta, phi,&
            &gYY)

        !GZZ
        CALL GetMetricComponent(&
            &M, Mr, NP,&
            &bufsize,&
            &timeofdata(7),&
            &iter, nchunks,&
            &6,&
            &gzzfilename,&
            &r, theta, phi,&
            &gZZ)

        !GXY
        CALL GetMetricComponent(&
            &M, Mr, NP,&
            &bufsize,&
            &timeofdata(8),&
            &iter, nchunks,&
            &7,&
            &gxyfilename,&
            &r, theta, phi,&
            &gXY)

        !GXZ
        CALL GetMetricComponent(&
            &M, Mr, NP,&
            &bufsize,&
            &timeofdata(9),&    
            &iter, nchunks,&
            &8,&
            &gxzfilename,&
            &r, theta, phi,&
            &gXZ)

        !GYZ
        CALL GetMetricComponent(&
            &M, Mr, NP,&
            &bufsize,&
            &timeofdata(10),&
            &iter, nchunks,&
            &9,&
            &gyzfilename,&
            &r, theta, phi,&
            &gYZ)


        DO i = 1, 9
            IF( timeofdata(i) .NE. timeofdata(10) ) STOP "Metric components time discrepancy."
        END DO

        time = timeofdata(10)

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
      END SUBROUTINE GetAllMetricComponents

