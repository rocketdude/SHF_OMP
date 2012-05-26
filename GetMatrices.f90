!--------------------------------------------------------!
!    Get Matrices subroutine                             !
!--------------------------------------------------------!

      SUBROUTINE GetMatrices(& 
& M, Mr, NP, Lmax, LWORK, IFLAG,&
& r, rho, theta, phi,&
& AF, AFinv, B, Dth, Dphi)

        !This subroutine initializes the following matrices
        !1. AF--matrix of Tn(r)Ylm(theta,phi) with {ijk} rows and {nlm} columns
        !   size(AF) = [4*np, np]
        !2. B--matrix of quadrature used to compute the values of the coefficients
        !   a_nlm(t) with {nlm} rows and {ijk} columns. size(B) = [np, 4*np]
        !3. Dth--matrices that are used to calculate the derivative
        !   of spherical harmonics wrt theta. size(Dth) = [4*np, np]

        !Dth is the matrix used to calculate S_{,theta}
        !by means of the recurrence relation:
        !Ylm_{,theta} = m*cot(theta)*Ylm + sqrt( (l-m)*(l+m+1) ) Yl,m+1 * exp(-i*phi)

        !4. Dphi--matrices that are used to calculate the derivatives
        !   of spherical harmonics wrt phi. size(Dphi) = [4*np, np]

        !Calculate matrix Dphi
        !Dphi is the matrix used to calculate S_{,phi} using:
        !Ylm_{,phi} = (im)*Ylm
        
        USE            omp_lib
        IMPLICIT       none

!--------------------------------------------------------!
!     Declare calling variables                          !
!--------------------------------------------------------!

        INTEGER*4, INTENT(in)  :: M, Mr, NP, Lmax, LWORK, IFLAG
        
        REAL*8, INTENT(in)   :: r(Mr+1)
        REAL*8, INTENT(in)   :: rho(Mr+1)
        REAL*8, INTENT(in)   :: theta(2*M)
        REAL*8, INTENT(in)   :: phi(2*M)

        COMPLEX*16, INTENT(out)  :: AF(4*NP, NP)
        COMPLEX*16, INTENT(out)  :: AFinv(NP, 4*NP)
        COMPLEX*16, INTENT(out)  :: B(NP, 4*NP)
        COMPLEX*16, INTENT(out)  :: Dth(4*NP, NP)
        COMPLEX*16, INTENT(out)  :: Dphi(4*NP, NP)

!--------------------------------------------------------!
!     Declare Locals                                     !
!--------------------------------------------------------!

        INTEGER*4        i, j, k
        INTEGER*4        n, l, ml
        INTEGER*4        crow      !Row counter
        INTEGER*4        ccol      !Row counter
        INTEGER*4      INFO

        CHARACTER*32   CTemp

        COMPLEX*16     f_nlm
        COMPLEX*16     g_nlm
        COMPLEX*16     f_nlmp1

        COMPLEX*16     AF2(4*NP, NP)
        COMPLEX*16     WORK(LWORK)
        COMPLEX*16     SIGMA(4*NP, NP)
        COMPLEX*16     U(4*NP, NP)
        COMPLEX*16     VT(NP, NP)

        REAL*8         RWORK(5*NP)
        REAL*8         S(NP)

!--------------------------------------------------------!
!      Main Subroutine                                   !
!--------------------------------------------------------!

        PRINT *, 'Initializing Matrices'

        !$OMP PARALLEL DO SHARED(rho, theta, phi, AF, B, Dth, Dphi)&
        !$OMP &PRIVATE(j, k, n, l, ml, crow, ccol, f_nlm, f_nlmp1, g_nlm)
        DO i = 1, (Mr+1)
           DO j = 1, 2*M
              DO k = 1, 2*M

                 crow = (i-1)*(2*M)*(2*M) + (j-1)*(2*M) + k
                 
                 DO n = 0, Mr
                    DO l = 0, Lmax
                       DO ml = -l, l

                          ccol = n*(Lmax+1)**2 + l**2 + (ml+l+1)

                          CALL EvaluateF(n, l, ml, rho(i), theta(j), phi(k), f_nlm)

                          IF( (ml+1) .GE. -l .AND. (ml+1) .LE. +l ) THEN
                             CALL EvaluateF(n, l, (ml+1), rho(i), theta(j), phi(k), f_nlmp1)
                          ELSE
                             f_nlmp1 = CMPLX(0.0D0, 0.0D0)
                          END IF

                          CALL EvaluateG(n, l, ml, rho(i), theta(j), phi(k), M, Mr, g_nlm)

                          !Populate Matrices
                          AF(crow, ccol) = f_nlm
                          B(ccol, crow) = g_nlm

                          Dth(crow, ccol) = &
                               &(DBLE(ml) / TAN( theta(j) ) *f_nlm)+&
                               &(SQRT( DBLE( (l-ml)*(l+ml+1) ) )*&
                               &EXP(CMPLX(0.0D0, -1.0D0*phi(k)))*f_nlmp1)

                          Dphi(crow, ccol) = CMPLX(0.0D0, 1.0D0)*DBLE(ml)*f_nlm

                       END DO
                    END DO
                 END DO

              END DO
           END DO
        END DO
        !$OMP END PARALLEL DO

        CALL OMP_SET_NUM_THREADS(1)
        
        CTemp = 'AFinv.dat'
        IF( IFLAG .LE. 0 ) THEN

           !Calculating the Moore-Penrose inverse of AF using SVD
           AF2 = AF

           !Compute the SVD of AF
           CALL ZGESVD( 'S', 'S', 4*NP, NP, AF2, 4*NP, S, U, 4*NP, VT, NP, WORK, LWORK, RWORK, INFO)

           !If LWORK == -1, query for the size of WORK matrix
           IF( LWORK .EQ. -1 ) THEN
              PRINT *, 'INFO=', INFO                      !If INFO = 0, then ZGELSVD ran properly
              PRINT *, 'Optimal size of WORK=', WORK(1)
              STOP
           END IF

           !Compute AFinv (PINV of AF), AFinv = VT**T * SIGMA * U**T. (**T means transpose)
           !$OMP PARALLEL DO
           DO j = 1, NP
              CALL ZSCAL( 4*NP, DCMPLX( 1.0 / S(j)), U(1,j), 1 )
           END DO
           !$OMP END PARALLEL DO

           CALL ZGEMM( 'C', 'C', NP, 4*NP, NP, 1.0D0, VT, NP, U, 4*NP,&
                &0.0D0, AFinv, NP)

           CALL OMP_SET_NUM_THREADS(16)
           
           IF( IFLAG .EQ. -1 ) THEN
              CALL Write2dC(NP, 4*NP, AFinv, CTemp)
           END IF

        ELSEIF (IFLAG .EQ. 1 ) THEN
           CALL Read2dC(NP, 4*NP, AFinv, CTemp)

        ELSE
           PRINT *, 'Specify the correct number for IFLAG:'
           PRINT *, 'IFLAG = 0 ==> Calculate AFinv (slow)'
           PRINT *, 'IFLAG = -1 ==> Calculate AFinv (slow) and write to a file ~ 1 GB'
           PRINT *, 'IFLAG = 1 ==> Read AFinv from a AFinv.dat'
           STOP
        END IF
           
        RETURN
      END SUBROUTINE GetMatrices
