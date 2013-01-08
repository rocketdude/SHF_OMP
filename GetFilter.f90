!--------------------------------------------------------!
!    GetFilter subroutine                                !
!--------------------------------------------------------!

      SUBROUTINE GetFilter(&
&Mr, filterP, eps,&
&Filter)
 
        USE             omp_lib
        IMPLICIT        none

!--------------------------------------------------------!
!     Declare calling variables                          !
!--------------------------------------------------------!

        INTEGER*4               :: Mr,filterP
        REAL*8                  :: eps
        REAL*8, INTENT(out)     :: Filter(Mr+1)

!--------------------------------------------------------!
!     Declare Locals                                     !
!--------------------------------------------------------!

        INTEGER*4         n
        
!--------------------------------------------------------!
!      Main Subroutine                                   !
!--------------------------------------------------------!

        !$OMP PARALLEL DO
        DO n = 0, Mr
            Filter(n+1) = EXP( LOG(eps)*(n/Mr)**(2*filterP) )
        END DO
        !$OMP END PARALLEL DO

        RETURN
    END SUBROUTINE GetFilter
