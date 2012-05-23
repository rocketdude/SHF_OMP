!----------------------------------------------------------!
!     GetMetric  subroutine                                !
!----------------------------------------------------------!

      SUBROUTINE GetMetric(&
&M, Mr, NP,&
&h10, h9,&
&nx10, ny10, nz10,&
&nx9, ny9, nz9,&
&Ox10, Oy10, Oz10,&
&Ox9, Oy9, Oz9,&
&Xmax10, Ymax10, Zmax10,&
&Xmax9, Ymax9, Zmax9,&
&r, theta, phi,&
&LM,&
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
        INTEGER*4, INTENT(in)  :: nx10, ny10, nz10
        INTEGER*4, INTENT(in)  :: nx9, ny9, nz9
        INTEGER*4, INTENT(in)  :: Ox10, Oy10, Oz10
        INTEGER*4, INTENT(in)  :: Ox9, Oy9, Oz9

        INTEGER*4, INTENT(in)  :: LM(64, 64)

        REAL*8, INTENT(in)   :: h10, h9
        REAL*8, INTENT(in)   :: Xmax10, Ymax10, Zmax10
        REAL*8, INTENT(in)   :: Xmax9, Ymax9, Zmax9

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



!----------------------------------------------------------!
!      Main                                                !
!----------------------------------------------------------!

        PRINT *, 'Getting Metric Data'

CALL  

        PRINT *, 'DONE!'

        RETURN
      END SUBROUTINE GetMetric
