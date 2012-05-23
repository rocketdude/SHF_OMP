!--------------------------------------------------------!
!    Initialize Filter subroutine                        !
!--------------------------------------------------------!

      SUBROUTINE InitializeFilter(& 
& M, Mr, NP, Lmax,&
& eps,&
& F)

        !This subroutine initializes the exponential filters
        !used for the Chebyshev polynomials

        USE            omp_lib
        IMPLICIT       none

!--------------------------------------------------------!
!     Declare calling variables                          !
!--------------------------------------------------------!

        INTEGER*4, INTENT(in)  :: M, Mr, NP, Lmax
        REAL*8, INTENT(in)     :: eps

        REAL*8, INTENT(out)  :: F(NP)

!--------------------------------------------------------!
!     Declare Locals                                     !
!--------------------------------------------------------!

        INTEGER*4        n, l, ml
        INTEGER*4        crow      !Row counter

!--------------------------------------------------------!
!      Main Subroutine                                   !
!--------------------------------------------------------!

        PRINT *, 'Initializing Filter'
                 
        !$OMP PARALLEL DO SHARED(F, eps, Mr, Lmax) PRIVATE(crow, l, ml, n)
        DO n = 0, Mr
           DO l = 0, Lmax
              DO ml = -l, l

                 crow = n*(Lmax+1)**2 + l**2 + (ml+l+1)
                 F(crow) = EXP( LOG(eps) * ( DBLE(n) / DBLE(Mr) )**8 )

              END DO
           END DO
        END DO
        !$OMP END PARALLEL DO

        PRINT *, 'DONE!'

        RETURN
      END SUBROUTINE InitializeFilter
