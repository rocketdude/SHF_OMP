!========================================================!
!    Cubic Spline Interpolation Subroutine               !
!========================================================!

!Cubic spline interpolation routines from Numerical Recipes
!in FORTRAN 77 by Press, Teukolsky, Vetterling, and Flannery
!To use this: call spline to compute y2, 
!             then use splint to interpolate

    SUBROUTINE ComputeSpline2ndDeriv(x,y,n,yp1,ypn,y2)  
    !This subroutine is used to compute y2 which is the second
    !derivatives of the inerpolating function at the tabulated
    !points x. Set yp1 and ypn to be 1e30 or larger to get
    !natural spline where the boundaries have zero second deriv.

      INTEGER*4     n,NMAX  
      PARAMETER     (NMAX=500)

      REAL*8        yp1,ypn,x(n),y(n),y2(n)  
      REAL*8        p,qn,sig,un,u(NMAX)
      INTEGER*4     i,k  

      if (yp1.gt..99D30) then  
        y2(1)=0.0D0
        u(1)=0.0D0
      else  
        y2(1)=-0.5D0 
        u(1)=(3.0D0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)  
      endif

      do i=2,n-1  
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))  
        p=sig*y2(i-1)+2.0D0 
        y2(i)=(sig-1.)/p  
        u(i)=(6.0D0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))&
             &/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p  
      end do

      if (ypn.gt..99e30) then  
        qn=0.0D0  
        un=0.0D0
      else  
        qn=0.5D0
        un=(3.0D0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))  
      endif

      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.0D0)  
      do k=n-1,1,-1  
        y2(k)=y2(k)*y2(k+1)+u(k)  
      end do

      return
    END SUBROUTINE

!=====================================================================!
    
    SUBROUTINE CubicSplineInterpolation(xa,ya,y2a,klo,n,x,y)
    !Given a value of x, the routine returns the cubic-spline
    !interpolated value y.
    !y2a is computed using subroutine ComputeSpline2ndDeriv

      INTEGER*4     n  
      INTEGER*4     khi,klo  
      REAL*8        x,y,xa(n),y2a(n),ya(n)  
      REAL*8        a,b,h  

      khi=klo+1

      h=xa(khi)-xa(klo)  
      if (h.eq.0.0D0) then 
        stop "***ERROR in SPLINT: x has to be monotonically increasing"
      end if

      a=(xa(khi)-x)/h  
      b=(x-xa(klo))/h  
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.

      return

    END SUBROUTINE


!=====================================================================!

    SUBROUTINE PolynomialInterpolation(xa,ya,n,x,y,dy)

      INTEGER*4     n,NMAX  
      PARAMETER (NMAX=10)

      INTEGER*4     i,m,ns  
      REAL*8        dy,x,y,xa(n),ya(n)  
      REAL*8        den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)

      ns=1  
      dif=abs(x-xa(1))

      do i=1,n  
        dift=abs(x-xa(i))  
        if (dift.lt.dif) then  
          ns=i  
          dif=dift  
        endif  
        c(i)=ya(i)  
        d(i)=ya(i)
      end do

      y=ya(ns)  
      ns=ns-1

      do m=1,n-1  
        do i=1,n-m  
          ho=xa(i)-x  
          hp=xa(i+m)-x  
          w=c(i+1)-d(i)  
          den=ho-hp  
          if(den.eq.0.) stop "ERROR in polint"
          den=w/den  
          d(i)=hp*den  
          c(i)=ho*den
        end do
        if (2*ns.lt.n-m)then  
          dy=c(ns+1)  
        else  
          dy=d(ns)  
          ns=ns-1  
        endif  
        y=y+dy  
      end do  
      
      return
    END SUBROUTINE
