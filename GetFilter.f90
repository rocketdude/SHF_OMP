!--------------------------------------------------------!
!    GetFilter subroutine                                !
!--------------------------------------------------------!

      SUBROUTINE GetFilter(&
&Lmax, filterP, eps,&
&Filter)
 
        USE             omp_lib
        IMPLICIT        none

!--------------------------------------------------------!
!     Declare calling variables                          !
!--------------------------------------------------------!

        INTEGER*4               :: Lmax,filterP
        REAL*8                  :: eps
        REAL*8, INTENT(out)     :: Filter(Lmax+1)

!--------------------------------------------------------!
!     Declare Locals                                     !
!--------------------------------------------------------!

        INTEGER*4         l
        
!--------------------------------------------------------!
!      Main Subroutine                                   !
!--------------------------------------------------------!

        !$OMP PARALLEL DO
        DO l = 0, Lmax
            !Filter(l+1) = EXP( LOG(eps)*(l/Lmax)**(2*filterP) )
            Filter(l+1) = 1.0D0
        END DO
        !$OMP END PARALLEL DO

        RETURN
    END SUBROUTINE GetFilter
