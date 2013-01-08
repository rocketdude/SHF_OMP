!--------------------------------------------------------!
!    ApplyFilter subroutine                              !
!--------------------------------------------------------!

      SUBROUTINE ApplyFilter(&
&Mr, Lmax,&
&Filter,&
&a)
 
        USE             omp_lib
        IMPLICIT        none

!--------------------------------------------------------!
!     Declare calling variables                          !
!--------------------------------------------------------!

        INTEGER*4               :: Mr,Lmax
        REAL*8                  :: Filter(Mr+1)
        COMPLEX*16              :: a(Mr+1,2,Lmax+1,Lmax+1)

!--------------------------------------------------------!
!     Declare Locals                                     !
!--------------------------------------------------------!

        INTEGER*4         n
        
!--------------------------------------------------------!
!      Main Subroutine                                   !
!--------------------------------------------------------!

        !$OMP PARALLEL DO
        DO n = 0, Mr
            a(n+1,:,:,:) = a(n+1,:,:,:) * Filter(n+1)
        END DO
        !$OMP END PARALLEL DO

        RETURN
    END SUBROUTINE ApplyFilter
