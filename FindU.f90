!========================================================!
!    FindU Subroutine                                    !
!========================================================!

     SUBROUTINE FindU(&
&M, Mr, NP,&
&r,&
&AF,a,&
&it, WriteSit,&
&U, Uave)

!This subroutine finds U(theta, phi) using cubic spline interpolation
!and then performing linear interpolation

        USE             omp_lib
        IMPLICIT        none

!-----------------------------------------------------------!
!       Declare calling variables                           !
!-----------------------------------------------------------!

       INTEGER*4            :: M, Mr, NP, it, WriteSit

       REAL*8               :: r(Mr+1)
       COMPLEX*16           :: AF(4*NP, NP)
       COMPLEX*16           :: a(NP)

       REAL*8, INTENT(out)  :: U(2*M, 2*M)
       REAL*8, INTENT(out)  :: Uave

!--------------------------------------------------------!
!     Declare Locals                                     !
!--------------------------------------------------------!

       CHARACTER*32      CTemp

       INTEGER*4, PARAMETER ::  SP=1000 !Number of points to be calculated on the spline
       INTEGER*4                i, j, k
       INTEGER*4                iindex, jkindex
       INTEGER*4                ilow, ilow2, flag

       COMPLEX*16             S(4*NP)

       REAL*8                 U_r(Mr+1)
       REAL*8                 U_r2(Mr+1) !Second derivatives of U_r at points along r

       REAL*8                 rr(SP)
       REAL*8                 UU(SP)
       REAL*8                 Temp
       
       REAL*8                 S100
       REAL*8                 deltar
       REAL*8                 Splint !Define the type for function Splint


!--------------------------------------------------------!
!      Main Subroutine                                   !
!--------------------------------------------------------!

       S100 = 100.0D0

       !Evaluate the eikonal data S
       S = MATMUL(AF, a)

       !Writing S into file (only at certain iterations)
       IF( MOD(it, WriteSit) .EQ. 0 ) CALL WriteS(4*NP, ABS(S), it)

       !$OMP PARALLEL DO PRIVATE(k, i, iindex, jkindex, U_r, U_r2, flag, ilow, ilow2, deltar, rr, UU)
       DO j = 1, 2*M
          DO k = 1, 2*M

             jkindex = (j-1)*(2*M) + k

             DO i = 1, (Mr+1)
                iindex = (i-1)*(4*M**2)
                U_r(i) = ABS( S(iindex + jkindex) )
             END DO

             !To find the radial distance of the light cone,
             !we fit cubic splines and then perform linear interpolation
             
             !Calculate the second derivatives and store it into U_r2

             CALL ComputeSpline2ndDeriv(r, U_r, Mr+1, 1.0D31, 1.0D31, U_r2) !Natural cubic spline

             !Find ilow and ilow+1 which are the indices that bound where S100 is located
             flag = 0
             DO i=(1+10), (Mr+1-10)
                IF( U_r(i) .GE. S100 ) THEN
                   ilow = i-1
                   flag = 1
                   EXIT
                END IF
             END DO

             IF( flag .EQ. 0) THEN           
                WRITE(*,*) 'Error in hunting'
                READ(*, '()')
             END IF

             !Calculate the equispaced points rr--between r(ilow) and r(ilow+1)
             !and, calculate the values of ABS(S) at those points
             deltar = (r(ilow+1) - r(ilow))/DBLE(SP-1)
             DO i = 1, SP
                rr(i) = r(ilow) + DBLE(i-1)*deltar
                CALL CubicSplineInterpolation(r, U_r, U_r2, ilow, Mr+1, rr(i), UU(i))

             END DO

             !Perform HUNT again, finding ilow and ilow+1 which are the indices that bound
             !where S100 is located
             flag = 0
             DO i=1, SP
                IF( UU(i) .GE. S100 ) THEN
                   ilow = i-1
                   flag = 1
                   EXIT
                END IF
             END DO

             IF( flag .EQ. 0) THEN           
                WRITE(*,*) 'Error in hunting'
                READ(*, '()')
             END IF

             !Find the value of r where S100 is located by performing linear interpolation
             U(j,k) = rr(ilow) + &
                  &( rr(ilow+1) - rr(ilow) )/( UU(ilow+1) - UU(ilow) ) *&
                  &(S100 - UU(ilow) )

          END DO
       END DO
       !$OMP END PARALLEL DO

       !Calculate the average U: Uave (only applies if there's obvious spherical symmetery)
       Uave = SUM( U )/ (4*M*M)

       RETURN
     END SUBROUTINE FindU
