!----------------------------------------------------------!
!     Get metric at current time subroutine                !
!----------------------------------------------------------!

    SUBROUTINE GetMetricAtCurrentTime(&
    &Nr, Nth, Nphi, TP,&
    &readdata, SFLAG,&
    &t, t_thresh, tdir,&
    &t_data, it_data,&    
    &nchunks,&
    &bufsize,&
    &rmaxX, rmaxY, rmaxZ,&
    &rminX, rminY, rminZ,&
    &rho, theta, phi,&
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

    INTEGER*4               Nr, Nth, Nphi, TP
    INTEGER*4               nchunks
    INTEGER*4               readdata
    INTEGER*4               SFLAG
    REAL*8                  t
    REAL*8                  tdir
    REAL*8                  t_thresh
    REAL*8                  t_data(TP)
    INTEGER*4               it_data(TP)

    INTEGER(HSIZE_T)        bufsize(3)

    REAL*8                  rmaxX, rmaxY, rmaxZ
    REAL*8                  rminX, rminY, rminZ
    REAL*8                  rho(Nr), theta(Nth), phi(Nphi)

    !Stored metric data
    REAL*8                  BgRR(TP,Nr,Nth,Nphi)
    REAL*8                  BgThTh(TP,Nr,Nth,Nphi)
    REAL*8                  BgPhiPhi(TP,Nr,Nth,Nphi)
    REAL*8                  BgRTh(TP,Nr,Nth,Nphi)
    REAL*8                  BgRPhi(TP,Nr,Nth,Nphi)
    REAL*8                  BgThPhi(TP,Nr,Nth,Nphi)

    REAL*8                  Balpha(TP,Nr,Nth,Nphi)
    REAL*8                  BbetaR(TP,Nr,Nth,Nphi)
    REAL*8                  BbetaTh(TP,Nr,Nth,Nphi)
    REAL*8                  BbetaPhi(TP,Nr,Nth,Nphi)

    !Metric at current time
    REAL*8                  gRR(Nr,Nth,Nphi)
    REAL*8                  gThTh(Nr,Nth,Nphi)
    REAL*8                  gPhiPhi(Nr,Nth,Nphi)
    REAL*8                  gRTh(Nr,Nth,Nphi)      
    REAL*8                  gRPhi(Nr,Nth,Nphi)
    REAL*8                  gThPhi(Nr,Nth,Nphi)

    REAL*8                  alpha(Nr,Nth,Nphi) 
    REAL*8                  betaR(Nr,Nth,Nphi)    
    REAL*8                  betaTh(Nr,Nth,Nphi)
    REAL*8                  betaPhi(Nr,Nth,Nphi)

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
            &Nr, Nth, Nphi,&
            &it_data(i), nchunks,&
            &bufsize,&
            &t_data(i),&
            &rmaxX, rmaxY, rmaxZ,&
            &rminX, rminY, rminZ,&
            &rho, theta, phi,&
            &Balpha(i,:,:,:),&
            &BbetaR(i,:,:,:), BbetaTh(i,:,:,:), BbetaPhi(i,:,:,:),&
            &BgRR(i,:,:,:), BgThTh(i,:,:,:), BgPhiPhi(i,:,:,:),&
            &BgRTh(i,:,:,:), BgRPhi(i,:,:,:), BgThPhi(i,:,:,:))

        END DO

        IF( tdir .LT. 0.0D0 ) THEN
        
            t_thresh = MINVAL(t_data) + &
            & (MAXVAL(t_data) - MINVAL(t_data))*DBLE(TP-2)/DBLE(2*(TP-1))

            IF(SFLAG .EQ. 0) t = MAXVAL(t_data)
    
        ELSEIF( tdir .GT. 0.0D0 ) THEN
        
            t_thresh = MAXVAL(t_data) - &
            & (MAXVAL(t_data) - MINVAL(t_data))*DBLE(TP-2)/DBLE(2*(TP-1))

            IF(SFLAG .EQ. 0) t = MINVAL(t_data)

        END IF

    ELSEIF (readdata .GT. 0 ) THEN
        
        i = readdata

        CALL GetAllMetricComponents(&
        &Nr, Nth, Nphi,&
        &it_data(i), nchunks,&
        &bufsize,&
        &t_data(i),&
        &rmaxX, rmaxY, rmaxZ,&
        &rminX, rminY, rminZ,&
        &rho, theta, phi,&
        &Balpha(i,:,:,:),&
        &BbetaR(i,:,:,:), BbetaTh(i,:,:,:), BbetaPhi(i,:,:,:),&
        &BgRR(i,:,:,:), BgThTh(i,:,:,:), BgPhiPhi(i,:,:,:),&
        &BgRTh(i,:,:,:), BgRPhi(i,:,:,:), BgThPhi(i,:,:,:))

        IF( tdir .LT. 0.0D0 ) THEN
        
            t_thresh = MINVAL(t_data) + &
            & (MAXVAL(t_data) - MINVAL(t_data))*DBLE(TP-2)/DBLE(2*(TP-1))
    
        ELSEIF( tdir .GT. 0.0D0 ) THEN
        
            t_thresh = MAXVAL(t_data) - &
            & (MAXVAL(t_data) - MINVAL(t_data))*DBLE(TP-2)/DBLE(2*(TP-1))

        END IF

    END IF

    !Check if there is any data with the same time
    DO j = 1, TP

        IF( ABS(t - t_data(j)) .LE. time_tolerance ) THEN
            sametime = j
        END IF

    END DO

    !If there's any data with the same time don't interpolate
    IF ( sametime .GT. 0 ) THEN

        k = sametime

        alpha = Balpha(k,:,:,:)
        betaR = BbetaR(k,:,:,:)
        betaTh = BbetaTh(k,:,:,:)
        betaPhi = BbetaPhi(k,:,:,:)
        gRR = BgRR(k,:,:,:)
        gThTh = BgThTh(k,:,:,:)
        gPhiPhi = BgPhiPhi(k,:,:,:)
        gRTh = BgRTh(k,:,:,:)
        gRPhi = BgRPhi(k,:,:,:)
        gThPhi = BgThPhi(k,:,:,:)

    ELSE

        DO i=1,Nr
          DO j=1,Nth
            DO k=1,Nphi

            CALL PolynomialInterpolation(&
            &t_data, Balpha(:,i,j,k),&
            &TP, t, alpha(i,j,k))

            CALL PolynomialInterpolation(&
            &t_data, BbetaR(:,i,j,k),&
            &TP, t, betaR(i,j,k))

            CALL PolynomialInterpolation(&
            &t_data, BbetaTh(:,i,j,k),&
            &TP, t, betaTh(i,j,k))

            CALL PolynomialInterpolation(&
            &t_data, BbetaPhi(:,i,j,k),&
            &TP, t, betaPhi(i,j,k)) 

            CALL PolynomialInterpolation(&
            &t_data, BgRR(:,i,j,k),&
            &TP, t, gRR(i,j,k))

            CALL PolynomialInterpolation(&
            &t_data, BgThTh(:,i,j,k),&
            &TP, t, gThTh(i,j,k))

            CALL PolynomialInterpolation(&
            &t_data, BgPhiPhi(:,i,j,k),&
            &TP, t, gPhiPhi(i,j,k))

            CALL PolynomialInterpolation(&
            &t_data, BgRTh(:,i,j,k),&
            &TP, t, gRTh(i,j,k))

            CALL PolynomialInterpolation(&
            &t_data, BgRPhi(:,i,j,k),&
            &TP, t, gRPhi(i,j,k))

            CALL PolynomialInterpolation(&
            &t_data, BgThPhi(:,i,j,k),&
            &TP, t, gThPhi(i,j,k))

            END DO
          END DO
        END DO

    END IF

    readdata = -1

    RETURN
    END SUBROUTINE GetMetricAtCurrentTime
