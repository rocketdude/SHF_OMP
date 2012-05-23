!---------------------------------------------------------!
!    FUNCTION Spherical Harmonic                          !
!---------------------------------------------------------!

      COMPLEX*16 FUNCTION SphHarmonicY(&
&l, ml,&
&theta, phi)

        IMPLICIT       none
        
        !This function calculates the value of Spherical Harmonic Y_lm
        !at a specific theta and phi
        !using the recurrence relation of the Legendre functions P_l^m
        !Refer to the Numerical Recipes in Fortran77 pg247 for more details

!---------------------------------------------------------!
!     Declare variables passed by function                !
!---------------------------------------------------------!

        INTEGER*4        l, ml
        REAL*8         theta, phi

        REAL*8        FactorialRatio

!---------------------------------------------------------!
!     Declare Locals                                      !
!---------------------------------------------------------!	

        INTEGER*4        i, ll   ! Counters for DO loops
        INTEGER*4        AbsM

        REAL*8         PI
        REAL*8         x       ! x = COS(theta)
        REAL*8         fact    ! odd factorial counter
        REAL*8         pmm     ! starting value P_m^m (note that P_(m-1)^m = 0 )
        REAL*8         pmmp1   ! P_(m+1)^m        
        REAL*8         pll     ! P_l^l
        REAL*8         somx2   ! somx2 = SQRT( 1 - x^2 )
        REAL*8         plgndr  ! plgndr is the value of P_l^m ( cos(theta) )
        REAL*8         epsilon ! sign factor in spherical harmonics, see Griffiths QM book

!---------------------------------------------------------!
!     Calculate SphHarmonicY                              !
!---------------------------------------------------------!
        
        PI = 3.141592653589793238462643383279502884197D0

        AbsM = ABS(ml)
        x = COS(theta)

        IF( AbsM .GT. l .or. l .LT. 0 ) THEN           
           WRITE(*,*) 'Bad arguments in SphHarmonicY'
           READ(*, '()')
        END IF

        pmm = 1.0D0

        IF( AbsM .GT. 0 ) THEN

           somx2 = SQRT( (1.0D0 - x) * (1.0D0 + x) )
           fact = 1.0D0

           DO i = 1, AbsM

              pmm = -1.0D0 * pmm * fact * somx2
              fact = fact + 2.0D0

           END DO

        END IF


        IF( l .EQ. AbsM ) THEN
           plgndr = pmm
        ELSE
           pmmp1 = x * DBLE(2*AbsM+1) * pmm
           
           IF( l .EQ. (AbsM+1) ) THEN
              plgndr = pmmp1
           ELSE
              DO ll = (AbsM+2), l
                 pll = (x * DBLE(2*ll - 1) * pmmp1 - DBLE(ll + AbsM - 1)*pmm )/DBLE(ll-AbsM)
                 pmm = pmmp1
                 pmmp1 = pll
              END DO

              plgndr = pll

           END IF

        END IF

        IF( ml .LT. 0 ) THEN
           epsilon = (-1.0D0)**AbsM
        ELSE
           epsilon = 1.0D0
        END IF

        SphHarmonicY = epsilon *&
             &SQRT( DBLE(2*l+1)/(4*PI)*FactorialRatio(l, AbsM) )*&
             &plgndr * CMPLX( COS(ml*phi), SIN(ml*phi) )

        RETURN
      END FUNCTION SphHarmonicY
