!----------------------------------------------------------!
!     GetMetric  subroutine                                !
!----------------------------------------------------------!

      SUBROUTINE GetAllMetricComponents(&
&Nr, Nth, Nphi,&
&iter, nchunks,&
&bufsize,&
&time,&
&rmax, rmin,&
&rho, theta, phi,&
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

        INTEGER*4               :: Nr, Nth, Nphi
        INTEGER*4               :: iter
        INTEGER*4               :: nchunks

        INTEGER(HSIZE_T)        :: bufsize(3)        

        REAL*8                  :: rmax, rmin
        REAL*8                  :: rho(Nr)
        REAL*8                  :: theta(Nth)
        REAL*8                  :: phi(Nphi)

        REAL*8, INTENT(out)     :: alpha(Nr,Nth,Nphi)
        REAL*8, INTENT(out)     :: betaR(Nr,Nth,Nphi)
        REAL*8, INTENT(out)     :: betaTh(Nr,Nth,Nphi)
        REAL*8, INTENT(out)     :: betaPhi(Nr,Nth,Nphi)

        REAL*8, INTENT(out)     :: gRR(Nr,Nth,Nphi)
        REAL*8, INTENT(out)     :: gThTh(Nr,Nth,Nphi)
        REAL*8, INTENT(out)     :: gPhiPhi(Nr,Nth,Nphi)

        REAL*8, INTENT(out)     :: gRTh(Nr,Nth,Nphi)
        REAL*8, INTENT(out)     :: gRPhi(Nr,Nth,Nphi)
        REAL*8, INTENT(out)     :: gThPhi(Nr,Nth,Nphi)

        REAL*8, INTENT(out)     :: time

!-----------------------------------------------------------!
!       Declare local variables                             !
!-----------------------------------------------------------!

        CHARACTER*32            format_string_Lapse(nchunks)
        CHARACTER*32            format_string_Shift(nchunks)
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

        REAL*8                  alphaXYZ(Nr,Nth,Nphi)
        REAL*8                  betaX(Nr,Nth,Nphi)
        REAL*8                  betaY(Nr,Nth,Nphi)
        REAL*8                  betaZ(Nr,Nth,Nphi)
        REAL*8                  gXX(Nr,Nth,Nphi)
        REAL*8                  gYY(Nr,Nth,Nphi)
        REAL*8                  gZZ(Nr,Nth,Nphi)
        REAL*8                  gXY(Nr,Nth,Nphi)
        REAL*8                  gXZ(Nr,Nth,Nphi)
        REAL*8                  gYZ(Nr,Nth,Nphi)

        REAL*8                  JMatrix(4,4)
        REAL*8                  gcart(4,4)
        REAL*8                  gsph(4,4)
        REAL*8                  dt_ratio
        REAL*8                  timeofdata(10)

        INTEGER*4               i, j, k
        INTEGER*4               error

!----------------------------------------------------------!
!      Main                                                !
!----------------------------------------------------------!

        PRINT *, 'Getting Metric Data of it=', iter

        !$OMP PARALLEL DO PRIVATE(j)
        DO i = 1,nchunks
            j = i - 1

            IF( j .LT. 10 ) THEN
                format_string_Lapse(i) = '(A19,I1,A3)'
                format_string_Shift(i) = '(A19,I1,A3)'
                format_string_SpatMetric(i) = '(A17,I1,A3)'
            ELSEIF( (j .GE. 10) .AND. (j .LT. 100) ) THEN
                format_string_Lapse(i) = '(A19,I2,A3)'
                format_string_Shift(i) = '(A19,I2,A3)'
                format_string_SpatMetric(i) = '(A17,I2,A3)'
            ELSEIF( (j .GE. 100) .AND. (j .LT. 1000) ) THEN
                format_string_Lapse(i) = '(A19,I3,A3)'
                format_string_Shift(i) = '(A19,I3,A3)'
                format_string_SpatMetric(i) = '(A17,I3,A3)'
            END IF

            IF( iter .LE. 8000 ) THEN
 
                WRITE(alphafilename(i), format_string_Lapse(i)) 't00_t50/alpha.file_', j, '.h5'
                WRITE(beta1filename(i), format_string_Shift(i)) 't00_t50/beta1.file_', j, '.h5'
                WRITE(beta2filename(i), format_string_Shift(i)) 't00_t50/beta2.file_', j, '.h5'
                WRITE(beta3filename(i), format_string_Shift(i)) 't00_t50/beta3.file_', j, '.h5'
                WRITE(gxxfilename(i), format_string_SpatMetric(i)) 't00_t50/gxx.file_', j, '.h5'
                WRITE(gyyfilename(i), format_string_SpatMetric(i)) 't00_t50/gyy.file_', j, '.h5'
                WRITE(gzzfilename(i), format_string_SpatMetric(i)) 't00_t50/gzz.file_', j, '.h5'
                WRITE(gxyfilename(i), format_string_SpatMetric(i)) 't00_t50/gxy.file_', j, '.h5'
                WRITE(gxzfilename(i), format_string_SpatMetric(i)) 't00_t50/gxz.file_', j, '.h5'
                WRITE(gyzfilename(i), format_string_SpatMetric(i)) 't00_t50/gyz.file_', j, '.h5'

            ELSE

                WRITE(alphafilename(i), format_string_Lapse(i)) 't50_t75/alpha.file_', j, '.h5'
                WRITE(beta1filename(i), format_string_Shift(i)) 't50_t75/beta1.file_', j, '.h5'
                WRITE(beta2filename(i), format_string_Shift(i)) 't50_t75/beta2.file_', j, '.h5'
                WRITE(beta3filename(i), format_string_Shift(i)) 't50_t75/beta3.file_', j, '.h5'
                WRITE(gxxfilename(i), format_string_SpatMetric(i)) 't50_t75/gxx.file_', j, '.h5'
                WRITE(gyyfilename(i), format_string_SpatMetric(i)) 't50_t75/gyy.file_', j, '.h5'
                WRITE(gzzfilename(i), format_string_SpatMetric(i)) 't50_t75/gzz.file_', j, '.h5'
                WRITE(gxyfilename(i), format_string_SpatMetric(i)) 't50_t75/gxy.file_', j, '.h5'
                WRITE(gxzfilename(i), format_string_SpatMetric(i)) 't50_t75/gxz.file_', j, '.h5'
                WRITE(gyzfilename(i), format_string_SpatMetric(i)) 't50_t75/gyz.file_', j, '.h5'

            END IF

        END DO
        !$OMP END PARALLEL DO
                
        !ALPHA
        CALL GetMetricComponent(&
            &Nr, Nth, Nphi,&
            &bufsize,&
            &timeofdata(1),&
            &iter, nchunks,&
            &0,&
            &alphafilename,&
            &rmax, rmin,&
            &rho, theta, phi,&
            &alphaXYZ)

        !BETA1
        CALL GetMetricComponent(&
            &Nr, Nth, Nphi,&
            &bufsize,&
            &timeofdata(2),&
            &iter, nchunks,&
            &1,&
            &beta1filename,&
            &rmax, rmin,&
            &rho, theta, phi,&
            &betaX)

        !BETA2
        CALL GetMetricComponent(&
            &Nr, Nth, Nphi,&
            &bufsize,&
            &timeofdata(3),&
            &iter, nchunks,&
            &2,&
            &beta2filename,&
            &rmax, rmin,&
            &rho, theta, phi,&
            &betaY)

        !BETA3
        CALL GetMetricComponent(&
            &Nr, Nth, Nphi,&
            &bufsize,&
            &timeofdata(4),&
            &iter, nchunks,&
            &3,&
            &beta3filename,&
            &rmax, rmin,&
            &rho, theta, phi,&
            &betaZ)

        !GXX
        CALL GetMetricComponent(&
            &Nr, Nth, Nphi,&
            &bufsize,&
            &timeofdata(5),&
            &iter, nchunks,&
            &4,&
            &gxxfilename,&
            &rmax, rmin,&
            &rho, theta, phi,&
            &gXX)

        !GYY
        CALL GetMetricComponent(&
            &Nr, Nth, Nphi,&
            &bufsize,&
            &timeofdata(6),&
            &iter, nchunks,&
            &5,&
            &gyyfilename,&
            &rmax, rmin,&
            &rho, theta, phi,&
            &gYY)

        !GZZ
        CALL GetMetricComponent(&
            &Nr, Nth, Nphi,&
            &bufsize,&
            &timeofdata(7),&
            &iter, nchunks,&
            &6,&
            &gzzfilename,&
            &rmax, rmin,&
            &rho, theta, phi,&
            &gZZ)

        !GXY
        CALL GetMetricComponent(&
            &Nr, Nth, Nphi,&
            &bufsize,&
            &timeofdata(8),&
            &iter, nchunks,&
            &7,&
            &gxyfilename,&
            &rmax, rmin,&
            &rho, theta, phi,&
            &gXY)

        !GXZ
        CALL GetMetricComponent(&
            &Nr, Nth, Nphi,&
            &bufsize,&
            &timeofdata(9),&    
            &iter, nchunks,&
            &8,&
            &gxzfilename,&
            &rmax, rmin,&
            &rho, theta, phi,&
            &gXZ)

        !GYZ
        CALL GetMetricComponent(&
            &Nr, Nth, Nphi,&
            &bufsize,&
            &timeofdata(10),&
            &iter, nchunks,&
            &9,&
            &gyzfilename,&
            &rmax, rmin,&
            &rho, theta, phi,&
            &gYZ)


        DO i = 1, 9
            IF( timeofdata(i) .NE. timeofdata(10) ) &
            & STOP "Metric components time discrepancy."
        END DO

        time = timeofdata(10)

        !Now we need to perform coordinate transformation 
        !from (t,x,y,z) to (t,rho,th,phi)
        !Both g_sph and g_cart are of down indices
 
        !$OMP PARALLEL DO PRIVATE(i, k, &
        !$OMP & JMatrix, gcart, gsph, error)
        DO j = 1, Nth
            DO k = 1, Nphi
                DO i = 1, Nr

                    !Build the Jacobian matrix JMatrix
                    CALL EvaluateJacobian(rmax,rmin,&
                                        &rho(i),theta(j),phi(k),&
                                        &JMatrix)
 
                    ! Build the matrix g_cart
                    CALL EvaluateMatrixofMetric(&
                         &alphaXYZ(i,j,k),&
                         &betaX(i,j,k), betaY(i,j,k), betaZ(i,j,k),&
                         &gXX(i,j,k), gYY(i,j,k), gZZ(i,j,k),&
                         &gXY(i,j,k), gXZ(i,j,k), gYZ(i,j,k),&
                         &gcart)

                    !  g_sph = (J^T) * g_cart * J
                    gsph = MATMUL(TRANSPOSE(JMatrix), MATMUL(gcart,JMatrix))

                    ! Get the 3-metric by inverting the 3-matric part of g_sph
                    CALL Invert3Metric(gsph(2,2), gsph(3,3), gsph(4,4),&
                                    &gsph(2,3), gsph(2,4), gsph(3,4),&
                                    &gRR(i,j,k), gThTh(i,j,k), gPhiPhi(i,j,k),&
                                    &gRTh(i,j,k), gRPhi(i,j,k), gThPhi(i,j,k),&
                                    &error)

                    IF( error .NE. 0 ) THEN
                        PRINT *, 'Determinant is zero, error in inverting metric'
                        PRINT *, 'i=',i,',j=',j,',k=',k
                        STOP
                    END IF

                    ! Get beta^r, beta^th, beta^phi
                    ! gsph(1,2) = g_01 = beta_r, gsph(1,3) = g_02 = beta_th, 
                    ! gsph(1,4) = g_03 = beta_phi
                    betaR(i,j,k) = gRR(i,j,k)*gsph(1,2) &
                                &+ gRTh(i,j,k)*gsph(1,3) &
                                &+ gRPhi(i,j,k)*gsph(1,4)
                    betaTh(i,j,k)= gRTh(i,j,k)*gsph(1,2) &
                                &+ gThTh(i,j,k)*gsph(1,3) &
                                &+ gThPhi(i,j,k)*gsph(1,4)
                    betaPhi(i,j,k)= gRPhi(i,j,k)*gsph(1,2) & 
                                &+ gThPhi(i,j,k)*gsph(1,3) & 
                                &+ gPhiPhi(i,j,k)*gsph(1,4)

                    ! Get alpha
                    ! alpha = sqrt( beta^i beta_i - g_00 )
                    ! beta^r = betaR, beta^th = betaTh, beta^phi = betaPhi, 
                    ! gsph(1,1) = g_00
                    alpha(i,j,k) = SQRT( gsph(1,2) * betaR(i,j,k) + &
                                    &gsph(1,3) * betaTh(i,j,k) + &
                                    &gsph(1,4) * betaPhi(i,j,k) &
                                & - gsph(1,1) )

                END DO
            END DO
        END DO
        !$OMP END PARALLEL DO

        PRINT *, 'DONE!'

        RETURN
      END SUBROUTINE GetAllMetricComponents

