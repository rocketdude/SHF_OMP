!--------------------------------------------------------!
!    Evaluate dSdr subroutine                            !
!--------------------------------------------------------!

      SUBROUTINE EvaluatedSdr(&
&M, Mr, NP, Lmax,&
&rmin, rmax,&
&a, AF,&
&dSdr)

        !This subroutine calculates S_{,r} at the collocation points using recursion.
        !The routine to calculate the derivative of the Chebyshev polynomials can be
        !found in Numerical Recipes in Fortran77.

        USE             omp_lib
        IMPLICIT        none

!--------------------------------------------------------!
!     Declare calling variables                          !
!--------------------------------------------------------!

        INTEGER*4, INTENT(in)  :: M, Mr, NP, Lmax

        REAL*8, INTENT(in)   :: rmin
        REAL*8, INTENT(in)   :: rmax

        COMPLEX*16, INTENT(in)  :: a(NP)
        COMPLEX*16, INTENT(in)  :: AF(4*NP, NP)

        COMPLEX*16, INTENT(out) :: dSdr(4*NP)

!--------------------------------------------------------!
!     Declare Locals                                     !
!--------------------------------------------------------!

        INTEGER*4         n
        INTEGER*4         nn
        INTEGER*4         nder
        INTEGER*4         cp

        REAL*8          const

        COMPLEX*16      D((Lmax+1)**2, Mr+1)
        COMPLEX*16      Dr((Lmax+1)**2, Mr+1)
        COMPLEX*16      ader(NP)
        COMPLEX*16      sum

!--------------------------------------------------------!
!      Main Subroutine                                   !
!--------------------------------------------------------!

        cp = (Lmax+1)**2
        const = 2.0D0/(rmax-rmin)
        
        !Rearrange vector a to matrix D
        !where rows correspond to {l,m} and columns correspond to {n}

        !$OMP PARALLEL DO
        DO n = 0, Mr
           D(:, n+1) = a( ((n*cp)+1) : ((n+1)*cp) )
        END DO
        !$OMP END PARALLEL DO

        !Find the derivatives by using recursion relations
        Dr(:, Mr+1) = (0.0D0, 0.0D0)
        Dr(:, Mr) = CMPLX(2.0D0*DBLE(Mr), 0.0D0)*D(:,Mr)

        DO nder = 1, (Mr-1)
           nn = Mr-nder+1
           Dr(:, nn-1) = Dr(:, nn+1) + CMPLX(2.0D0*DBLE(nn-1),0.0D0)*D(:,nn)
        END DO

        Dr = const * Dr

        DO n = 0, Mr
           ader( ((n*cp)+1) : ((n+1)*cp) ) = Dr(:, n+1)
        END DO

        dSdr = MATMUL(AF, ader)
        RETURN
      END SUBROUTINE EvaluatedSdr
