!----------------------------------------------------------!
!     Get metric at current time subroutine                !
!----------------------------------------------------------!

    SUBROUTINE GetMetricAtCurrentTime(&
    &NP, TP,&
    &t, t_data,&
    &alpha, betaR, betaTh, betaPhi,&
    &gRR, gThTh, gPhiPhi,&
    &gRTh, gRPhi, gThPhi,&
    &Balpha, BbetaR, BbetaTh, BbetaPhi,&
    &BgRR, BgThTh, BgPhiPhi,&
    &BgRTh, BgRPhi, BgThPhi)



    USE omp_lib
    IMPLICIT  none

!-----------------------------------------------------------!
!       Declare calling variables                           !
!-----------------------------------------------------------!

    INTEGER*4     NP
    REAL*8        t, tA, tB

    !Metric data at time tA --> will be changed to metric
    !data at time t and tA will be equal to t
    REAL*8        gRR(4*NP)
    REAL*8        gThTh(4*NP)
    REAL*8        gPhiPhi(4*NP)
    REAL*8        gRTh(4*NP)      
    REAL*8        gRPhi(4*NP)
    REAL*8        gThPhi(4*NP)

    REAL*8        alpha(4*NP) 
    REAL*8        betaR(4*NP)    
    REAL*8        betaTh(4*NP)
    REAL*8        betaPhi(4*NP)

    !Metric data at time tB
    REAL*8        BgRR(4*NP)
    REAL*8        BgThTh(4*NP)
    REAL*8        BgPhiPhi(4*NP)
    REAL*8        BgRTh(4*NP)
    REAL*8        BgRPhi(4*NP)
    REAL*8        BgThPhi(4*NP)

    REAL*8        Balpha(4*NP)
    REAL*8        BbetaR(4*NP)
    REAL*8        BbetaTh(4*NP)
    REAL*8        BbetaPhi(4*NP)

!-----------------------------------------------------------!
!       Declare local variables                             !
!-----------------------------------------------------------!

    REAL*8        IntRatio      !Interpolation ratio
    
!-----------------------------------------------------------!
!      Main                                                 !
!-----------------------------------------------------------!

    IF( tB .EQ. tA ) THEN
        STOP '***Error tB = tA, cannot interpolate***'
    ELSE
        IntRatio = (t-tA)/(tB-tA)
        
        alpha = alpha + (Balpha - alpha) * IntRatio
        betaR = betaR + (BbetaR - betaR) * IntRatio
        betaTh = betaTh + (BbetaTh - betaTh) * IntRatio
        betaPhi = betaPhi + (BbetaPhi - betaPhi) * IntRatio
        gRR = gRR + (BgRR - gRR) * IntRatio
        gThTh = gThTh + (BgThTh - gThTh) * IntRatio
        gPhiPhi = gPhiPhi + (BgPhiPhi - gPhiPhi) * IntRatio
        gRTh = gRTh + (BgRTh - gRTh) * IntRatio
        gRPhi = gRPhi + (BgRPhi - gRPhi) * IntRatio
        gThPhi = gThPhi + (BgThPhi - gThPhi) * IntRatio

        tA = t
    END IF
    
    RETURN
    END SUBROUTINE TemporalMetricInterpolation
