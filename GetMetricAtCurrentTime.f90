!----------------------------------------------------------!
!     Get metric at current time subroutine                !
!----------------------------------------------------------!

    SUBROUTINE GetMetricAtCurrentTime(&
    &M, Mr, NP, TP,&
    &readdata,&
    &t, t_thresh, tdir,&
    &t_data, it_data,&    
    &nchunks,&
    &bufsize,&
    &r, theta, phi,&
    &Balpha, BbetaR, BbetaTh, BbetaPhi,&
    &BgRR, BgThTh, BgPhiPhi,&
    &BgRTh, BgRPhi, BgThPhi,&
    &alpha, betaR, betaTh, betaPhi,&
    &gRR, gThTh, gPhiPhi,&
    &gRTh, gRPhi, gThPhi)

    !Get the metric at current time by polynomial interpolation.

    USE omp_lib
    USE HDF5
    IMPLICIT  none

!-----------------------------------------------------------!
!       Declare calling variables                           !
!-----------------------------------------------------------!

    INTEGER*4               M, Mr, NP, TP
    INTEGER*4               nchunks
    INTEGER*4               readdata
    REAL*8                  t
    REAL*8                  t_thresh
    REAL*8                  t_data(TP)
    INTEGER*4               it_data(TP)

    INTEGER(HSIZE_T)        bufsize(3)

    REAL*8                  r(Mr+1), theta(2*M), phi(2*M)

    !Stored metric data
    REAL*8                  BgRR(TP, 4*NP)
    REAL*8                  BgThTh(TP, 4*NP)
    REAL*8                  BgPhiPhi(TP, 4*NP)
    REAL*8                  BgRTh(TP, 4*NP)
    REAL*8                  BgRPhi(TP, 4*NP)
    REAL*8                  BgThPhi(TP, 4*NP)

    REAL*8                  Balpha(TP, 4*NP)
    REAL*8                  BbetaR(TP, 4*NP)
    REAL*8                  BbetaTh(TP, 4*NP)
    REAL*8                  BbetaPhi(TP, 4*NP)

    !Metric at current time
    REAL*8                  gRR(4*NP)
    REAL*8                  gThTh(4*NP)
    REAL*8                  gPhiPhi(4*NP)
    REAL*8                  gRTh(4*NP)      
    REAL*8                  gRPhi(4*NP)
    REAL*8                  gThPhi(4*NP)

    REAL*8                  alpha(4*NP) 
    REAL*8                  betaR(4*NP)    
    REAL*8                  betaTh(4*NP)
    REAL*8                  betaPhi(4*NP)

!-----------------------------------------------------------!
!       Declare local variables                             !
!-----------------------------------------------------------!

    INTEGER*4               i, j, k
    INTEGER*4               sametime
    REAL*8                  time_tolerance

!-----------------------------------------------------------!
!      Main                                                 !
!-----------------------------------------------------------!

    time_tolerance = 1.0D-9
    sametime = -1

    !Read metric data when necessary
    IF (readdata .EQ. 0 ) THEN

        DO i = 1, TP
        
            CALL GetAllMetricComponents(&
            &M, Mr, NP,&
            &it_data(i), nchunks,&
            &bufsize,&
            &t_data(i),&
            &r, theta, phi,&
            &Balpha(i,:),&
            &BbetaR(i,:), BbetaTh(i,:), BbetaPhi(i,:),&
            &BgRR(i,:), BgThTh(i,:), BgPhiPhi(i,:),&
            &BgRTh(i,:), BgRPhi(i,:), BgThPhi(i,:))

        END DO

        IF( tdir .LT. 0.0D0 ) THEN
        
            t = MAXVAL(t_data)
            t_thresh = MINVAL(t_data) + &
                    & (MAXVAL(t_data) - MINVAL(t_data))*DBLE(TP-2)/DBLE(2*(TP-1))
    
        ELSEIF( tdir .GT. 0.0D0 ) THEN
        
            t = MINVAL(t_data)
            t_thresh = MAXVAL(t_data) - &
                    & (MAXVAL(t_data) - MINVAL(t_data))*DBLE(TP-2)/DBLE(2*(TP-1))

        END IF

    ELSEIF (readdata .GT. 0 ) THEN
        
        i = readdata

        CALL GetAllMetricComponents(&
        &M, Mr, NP,&
        &it_data(i), nchunks,&
        &bufsize,&
        &t_data(i),&
        &r, theta, phi,&
        &Balpha(i,:),&
        &BbetaR(i,:), BbetaTh(i,:), BbetaPhi(i,:),&
        &BgRR(i,:), BgThTh(i,:), BgPhiPhi(i,:),&
        &BgRTh(i,:), BgRPhi(i,:), BgThPhi(i,:))

        IF( tdir .LT. 0.0D0 ) THEN
        
            t_thresh = MINVAL(t_data) + &
                    & (MAXVAL(t_data) - MINVAL(t_data))*DBLE(TP-2)/DBLE(2*(TP-1))
    
        ELSEIF( tdir .GT. 0.0D0 ) THEN
        
            t_thresh = MAXVAL(t_data) - &
                    & (MAXVAL(t_data) - MINVAL(t_data))*DBLE(TP-2)/DBLE(2*(TP-1))

        END IF

    END IF

    !Check if there is any data with the same time
    DO j = 1, NP

        IF( ABS(t - t_data(j)) .LE. time_tolerance ) THEN
            sametime = j
        END IF

    END DO

    !If there's any data with the same time don't interpolate
    IF ( sametime .GT. 0 ) THEN

        k = sametime

        alpha = Balpha(k,:)
        betaR = BbetaR(k,:)
        betaTh = BbetaTh(k,:)
        betaPhi = BbetaPhi(k,:)
        gRR = BgRR(k,:)
        gThTh = BgThTh(k,:)
        gPhiPhi = BgPhiPhi(k,:)
        gRTh = BgRTh(k,:)
        gRPhi = BgRPhi(k,:)
        gThPhi = BgThPhi(k,:)

    ELSE

        DO k = 1, NP

            CALL PolynomialInterpolation(&
            &t_data, Balpha(:,k),&
            &TP, t, alpha(k))

            CALL PolynomialInterpolation(&
            &t_data, BbetaR(:,k),&
            &TP, t, betaR(k))

            CALL PolynomialInterpolation(&
            &t_data, BbetaTh(:,k),&
            &TP, t, betaTh(k))

            CALL PolynomialInterpolation(&
            &t_data, BbetaPhi(:,k),&
            &TP, t, betaPhi(k)) 

            CALL PolynomialInterpolation(&
            &t_data, BgRR(:,k),&
            &TP, t, gRR(k))

            CALL PolynomialInterpolation(&
            &t_data, BgThTh(:,k),&
            &TP, t, gThTh(k))

            CALL PolynomialInterpolation(&
            &t_data, BgPhiPhi(:,k),&
            &TP, t, gPhiPhi(k))

            CALL PolynomialInterpolation(&
            &t_data, BgRTh(:,k),&
            &TP, t, gRTh(k))

            CALL PolynomialInterpolation(&
            &t_data, BgRPhi(:,k),&
            &TP, t, gRPhi(k))

            CALL PolynomialInterpolation(&
            &t_data, BgThPhi(:,k),&
            &TP, t, gThPhi(k))

        END DO

    END IF

    readdata = -1

    RETURN
    END SUBROUTINE GetMetricAtCurrentTime
