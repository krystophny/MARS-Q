      subroutine spline2d(fun,x,y,nx,ny,kx,coef)
c-----------------------------------------------------------------------
c     setup routine for bicubic spline of fun[x,y]
c     output of this routine is nx*ny spline coefficients 
c     stored in coef(kx,ny)
c-----------------------------------------------------------------------
      implicit none
      real fun(*),x(*),y(*),coef(*)
      integer kx,nx,ny
      integer ind,j
      do j=1,ny
         ind=(j-1)*kx+1
         call spline(x,fun(ind),nx,-1.e30,-1.e30,coef(ind))
      end do
      return
      end
c=======================================================================
      subroutine spline2dt(fun_new,x_new,y_new,kx_new,nx_new,ny_new,
     $     fun_old,x_old,y_old,kx_old,nx_old,ny_old,coef)
c-----------------------------------------------------------------------
c     evaluates bicubic spline: determines fun_new[x_new,y_new]
c     after spline2d has been called to determine coef
c     in the calling program fun_new,run_old, and coef are 2d arrays:
c     fun_old(kx_old,ky_old),fun_new(kx_new,ky_new),coef(kx_old,ky_old)
c-----------------------------------------------------------------------
      implicit none
      real fun_new(*),x_new(*),y_new(*)
      real fun_old(*),x_old(*),y_old(*)
      real coef(*)
      integer kx_new,nx_new,ny_new
      integer kx_old,nx_old,ny_old
      integer i,j,ky_old,ind
      parameter(ky_old=6000)
      real ftemp(ky_old),ctemp(ky_old)
      if(ny_old.gt.ky_old) then
         write(6,'("dimensioning problem in spline2dt")')
         stop
      end if
      do i=1,nx_new
         do j=1,ny_old
            ind=(j-1)*kx_old+1
            call splint(x_old,fun_old(ind),coef(ind),nx_old,x_new(i),
     $           ftemp(j))
         end do
         call  spline(y_old,ftemp,ny_old,-1.e30,-1.e30,ctemp)
         do j=1,ny_new
            call splint(y_old,ftemp,ctemp,ny_old,y_new(j),
     $           fun_new((j-1)*kx_new+i))
         end do
      end do
      return
      end
c=======================================================================
      subroutine spline1d(ynew,xnew,nnew,yold,xold,nold,y2old)
c-----------------------------------------------------------------------
c     use 1d cubic spline on yold[xold] to produce ynew[xnew]
c     y2old(1:nold) is a work array
c     ynew(1:nnew) is the output
c-----------------------------------------------------------------------
      implicit none
      real ynew(*),yold(*),xnew(*),xold(*)
      real y2old(*)
      real yp1,ypn
      integer nnew,nold,i
      yp1=-1.e30
      ypn=-1.e30
      call spline(xold,yold,nold,yp1,ypn,y2old)
      do i=1,nnew
         call splint(xold,yold,y2old,nold,xnew(i),ynew(i))
c         write (*,'("i,xnew,ynew",i4,1p2e12.4)')i,xnew(i),ynew(i)
      end do
      return
      end
c=======================================================================
      subroutine spline1dr(ypnew,xnew,y0,x0,yold,xold,nold,y2old)
c-----------------------------------------------------------------------
c     use 1d cubic spline on yold[xold] to find a root y(xnew)=y0
c     with the initial guess x0, and the first derivative ypnew
c     at the root xnew.
c     y2old(1:nold) is a work array
c     ynew and xnew are the output
c-----------------------------------------------------------------------
      implicit none
      real yold(*),xold(*)
      real y2old(*)
      real yp1,ypn,x0,y0,xnew,ypnew
      integer nold
      yp1=-1.e30
      ypn=-1.e30
      call spline(xold,yold,nold,yp1,ypn,y2old)
      call rsplint(xold,yold,y2old,nold,x0,y0,xnew,ypnew)
      return
      end
c=======================================================================
      subroutine spline1dr_safe(ypnew,xnew,y0,x0,yold,xold,nold,
     >     y2old,ibracket,istat)
c-----------------------------------------------------------------------
c     Safeguarded cubic-spline root solve on one prescribed bracket.
c     The bracket is never left.  Status 1 is a strict interior root;
c     status 2/3 is an exact lower/upper mesh node and is rejected by
c     KLAM0_SAFE; status 0 is not a sign-changing bracket; status 4 is
c     an iteration failure.  x0 is only an initial point inside the
c     bracket; bisection remains the governing update.
c-----------------------------------------------------------------------
      implicit none
      real*8 yold(*),xold(*),y2old(*)
      real*8 ypnew,xnew,y0,x0
      real*8 xlo,xhi,xmid,flo,fhi,fmid,ymid,ypmid
      real*8 ftol,xtol,scale
      integer nold,ibracket,istat,iter

      istat=0
      ypnew=0.
      xnew=x0
      if (nold.lt.2) return
      if (ibracket.lt.1.or.ibracket.ge.nold) return

      call spline(xold,yold,nold,-1.e30,-1.e30,y2old)
      xlo=xold(ibracket)
      xhi=xold(ibracket+1)
      if (xhi.le.xlo) return
      scale=max(1d0,abs(y0),abs(yold(ibracket)),
     >          abs(yold(ibracket+1)))
      ftol=1d-12*scale
      xtol=1d-12*max(1d0,abs(xlo),abs(xhi))
      call spline1dr_eval(xold,yold,y2old,ibracket,xlo,flo,ypnew)
      flo=flo-y0
      call spline1dr_eval(xold,yold,y2old,ibracket,xhi,fhi,ypnew)
      fhi=fhi-y0
      if (abs(flo).le.ftol) then
         xnew=xlo
         istat=2
         call spline1dr_eval(xold,yold,y2old,ibracket,xnew,
     >                       fmid,ypnew)
         return
      endif
      if (abs(fhi).le.ftol) then
         xnew=xhi
         istat=3
         call spline1dr_eval(xold,yold,y2old,ibracket,xnew,
     >                       fmid,ypnew)
         return
      endif
      if (flo*fhi.ge.0.) return

      xmid=0.5*(xlo+xhi)
      if (x0.gt.xlo.and.x0.lt.xhi) xmid=x0
      do iter=1,128
         if (xmid.le.xlo.or.xmid.ge.xhi) xmid=0.5*(xlo+xhi)
         call spline1dr_eval(xold,yold,y2old,ibracket,xmid,
     >                       ymid,ypmid)
         fmid=ymid-y0
         if (abs(fmid).le.ftol.or.abs(xhi-xlo).le.xtol) then
            xnew=xmid
            ypnew=ypmid
            istat=1
            return
         endif
         if (flo*fmid.lt.0.) then
            xhi=xmid
            fhi=fmid
         else
            xlo=xmid
            flo=fmid
         endif
         xmid=0.5*(xlo+xhi)
      enddo
      xnew=xmid
      ypnew=ypmid
      istat=4
      return
      end
c=======================================================================
      subroutine spline1dr_eval(xa,ya,y2a,k,x,y,yp)
c-----------------------------------------------------------------------
c     Evaluate a cubic spline and its derivative on interval k.
c-----------------------------------------------------------------------
      implicit none
      real*8 xa(*),ya(*),y2a(*),x,y,yp
      real*8 h,a,b
      integer k

      h=xa(k+1)-xa(k)
      if (h.eq.0.) then
         y=0.
         yp=0.
         return
      endif
      a=(xa(k+1)-x)/h
      b=(x-xa(k))/h
      y=a*ya(k)+b*ya(k+1)+
     >  ((a**3-a)*y2a(k)+(b**3-b)*y2a(k+1))*(h**2)/6.
      yp=(ya(k+1)-ya(k))/h-(3.*a**2-1.)/6.*h*y2a(k)
     >   +(3.*b**2-1.)/6.*h*y2a(k+1)
      return
      end
c=======================================================================
      subroutine spline(x,y,n,yp1,ypn,y2)
c-----------------------------------------------------------------------
c     spline routine based upon numerical recipes
c     this is the setup routine which needs to be called only once
c     splines y as a function of x--both arrays have n elements
c     yp1 and ypn are boundary conditions on the spline
c     yp1=y'[x] at x=x[1]
c     if yp1>=1.e30 then y''[x]=0 at x=x[1] is used
c     if yp1<=-1.e30 then y'[x[1]] is calculated from first four points
c     ypn=y'[x] at x=x[n]
c     if ypn>=1.e30 then y''[x]=0 at x=x[n] is used
c     if ypn<=-1.e30 then y'[x[n]] is calculated from last four points
c     y2[1:n] is calculated array of the second derivatives of the
c     interpolating function at the x[i]
c-----------------------------------------------------------------------
      parameter (nmax=60000)
      real*8::x(n),y(n),y2(n),u(nmax)
      real*8::yp1t,ypnt,yp1,ypn
      if (n .gt. nmax)then
         stop 'spline;  dimensional error '
      endif
      a = .99e30
      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else if(yp1.lt.-.99e30) then
         yp1t=(3*x(1)**2+x(2)*x(3)+
     $        x(2)*x(4)+x(3)*x(4)-2*x(1)*(x(2)+x(3)+x(4)))*
     $ y(1)/((x(1)-x(2))*(x(1)-x(3))*(x(1)-x(4)))
         yp1t=yp1t+
     $ (-x(1)**2+x(1)*x(3)+x(1)*x(4)-x(3)*x(4))*y(2)/
     $ ((x(1)-x(2))*(x(2)-x(3))*(x(2)-x(4)))+
     $ (x(1)**2-x(1)*x(2)-x(1)*x(4)+x(2)*x(4))*y(3)/
     $ ((x(1)-x(3))*(x(2)-x(3))*(x(3)-x(4)))+
     $ (-x(1)**2+x(1)*x(2)+x(1)*x(3)-x(2)*x(3))*y(4)/
     $ ((x(1)-x(4))*(x(2)-x(4))*(x(3)-x(4)))
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1t)
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do  i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     *      /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      end do
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else if(ypn.lt.-.99e30) then
         ypnt=(-(x(-2+n)*x(-1+n))+x(-2+n)*x(n)+x(-1+n)*x(n)-x(n)**2)*
     $  y(-3+n)/
     $  ((-x(-3+n)+x(-2+n))*(-x(-3+n)+x(-1+n))*
     $  (-x(-3+n)+x(n)))
         ypnt=ypnt+
     $  (x(-3+n)*x(-1+n)-x(-3+n)*x(n)-x(-1+n)*x(n)+x(n)**2)*
     $  y(-2+n)/
     $  ((-x(-3+n)+x(-2+n))*(-x(-2+n)+x(-1+n))*
     $  (-x(-2+n)+x(n)))
         ypnt=ypnt+
     $  (-(x(-3+n)*x(-2+n))+x(-3+n)*x(n)+x(-2+n)*x(n)-x(n)**2)*
     $  y(-1+n)/
     $  ((-x(-3+n)+x(-1+n))*(-x(-2+n)+x(-1+n))*
     $  (-x(-1+n)+x(n)))
         ypnt=ypnt+
     $        (x(-3+n)*x(-2+n)+x(-3+n)*x(-1+n)+x(-2+n)*x(-1+n)-
     $        2*(x(-3+n)+x(-2+n)+x(-1+n))*x(n)+3*x(n)**2)*y(n)/
     $        ((-x(-3+n)+x(n))*(-x(-2+n)+x(n))*(-x(-1+n)+x(n)))
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypnt-(y(n)-y(n-1))/(x(n)-x(n-1)))
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do  k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
      end do
      return
      end
c=======================================================================
      subroutine splint(xa,ya,y2a,n,x,y)
c-----------------------------------------------------------------------
c     cubic spline evaluator--spline must be called first to evaluate
c     y2a
c     ya is a function of xa--both are arrays of length n
c     ya2[1:n] contains spline coefficients calculated in spline
c     x is the argument of y[x] where y is to be evaluated
c     y=y[x] is the returned value
c-----------------------------------------------------------------------
      dimension xa(n),ya(n),y2a(n)
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) then
         stop 'bad xa input.'
      endif
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+
     *      ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      return
      end
c=======================================================================
      subroutine zspline(xa,ya,y2a,n,zg,ng,za)
c-----------------------------------------------------------------------
c     zspline integrates the cubic spline of ya[xa]
c     it assumes that spline has already been called to evaluate y2a
c     xa[1:n], ya[1:n], y2a[1:n]
c     In Mathematica notation z[i]=zg+Integrate[y[x],{x,x[ng],x[i]}]
c     where both zg and ng are input quantities,
c     za is calculated here. 
c     za can then be used in zsplint to determine z at a 
c     specific x
c-----------------------------------------------------------------------
      implicit none
c
      real xa(*),ya(*),y2a(*),zg,za(*)
      integer n,ng
c
      integer j
      real const
c
      if (ng .lt. 0 .or. ng .gt. n)stop 'zspline: wrong ng'
c
      za(1)=0.
      do j=2,n
         za(j)=za(j-1)+0.5*(xa(j)-xa(j-1))*(ya(j)+ya(j-1))-
     >        (xa(j)-xa(j-1))**3*(y2a(j)+y2a(j-1))/24.0
      enddo
c
      const=zg-za(ng)
      do j=1,n
         za(j)=za(j)+const
      enddo
c
      return
      end
c=======================================================================
      subroutine zsplint(xa,ya,y2a,za,n,x,y,yp,z)
c-----------------------------------------------------------------------
c     evaluate cubic spline to determine function (y),
c     derivative (yp) and integral (z) at location x.
c     first spline must be called to obtain y2a and
c     zspline must be called to obtain za
c-----------------------------------------------------------------------
      dimension xa(n),ya(n),y2a(n),za(n)
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) then
         stop 'bad xa input.'
      endif
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+
     *     ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      yp=(ya(khi)-ya(klo))/h-(3.*a**2-1.)/6.*h*y2a(klo)
     >     +(3.*b**2-1.)/6.*h*y2a(khi)
      z=za(klo)+0.5*h*ya(klo)*(1.-a**2)+0.5*ya(khi)*h*b**2+
     >     y2a(klo)*h**3*(2.*a**2-a**4-1.)/24.-
     >     y2a(khi)*h**3*(2.*b**2-b**4)/24.
      return
      end
c=======================================================================
      subroutine rsplint(xa,ya,y2a,n,x0,y0,x,yp)
c-----------------------------------------------------------------------
c     find root of the cubic spline y(x)=y0 and the first derivative 
c     (yp) at the root, giving an initial guess x0 being a scalar. 
c     the initial guess x0 must be such that both x0 and x lie in the 
c     same inteval [xa(klo),xa(khi)], otherwise the procedure fails.
c     Newton iteration method is used to search the root, thus 
c     iteration may fail if y(x) is not monotonic at the interval.
c     first spline must be called to obtain y2a
c-----------------------------------------------------------------------
      dimension xa(n),ya(n),y2a(n)
      eps = 1.0e-12
      xo=x0
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.xo)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) then
         stop 'bad xa input.'
      endif
      y=y0+eps*10
      k=0
 2    if (abs(y-y0).gt.eps.and.k.lt.100) then
        a=(xa(khi)-xo)/h
        b=(xo-xa(klo))/h
        y=a*ya(klo)+b*ya(khi)+
     *     ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
        yp=(ya(khi)-ya(klo))/h-(3.*a**2-1.)/6.*h*y2a(klo)
     >     +(3.*b**2-1.)/6.*h*y2a(khi)
        x=(y0-y)/yp+xo
        k=k+1
        xo=x
        goto 2
      endif
      if (k.eq.100) then
        write(*,*) 'check rsplint: x0,y0=',x0,y0
        do k=1,n
           write(*,*) k,xa(k),ya(k)
        enddo
        write(*,*) 'newton iteration...'
        xo=x0
        y=y0+eps*10
        do k=1,100
        a=(xa(khi)-xo)/h
        b=(xo-xa(klo))/h
        y=a*ya(klo)+b*ya(khi)+
     *     ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
        yp=(ya(khi)-ya(klo))/h-(3.*a**2-1.)/6.*h*y2a(klo)
     >     +(3.*b**2-1.)/6.*h*y2a(khi)
        x=(y0-y)/yp+xo
        write(*,*) xo,y,x
        xo=x
        enddo
        write(*,*) 'rsplint:bad x0 input.'
      endif
      return
      end

c=======================================================================
      function zspline1d(ya,xa,n,y2,yy)
c-----------------------------------------------------------------------
c     use 1d cubic spline on ya[xa] to do integration
c     y2(1:n) is a work array
c-----------------------------------------------------------------------
      implicit none
      real    xa(*),y2(*),yy(*)
      real    yp1,ypn,zr,zi
      complex ya(*)
      complex zspline1d 
      integer n,j
      yp1=-1.e30
      ypn=-1.e30
      do j=1,n
         yy(j) = real(ya(j))
      enddo
      call spline(xa,yy,n,yp1,ypn,y2)
      zr=0.
      do j=2,n
         zr=zr+0.5*(xa(j)-xa(j-1))*(yy(j)+yy(j-1))-
     >      (xa(j)-xa(j-1))**3*(y2(j)+y2(j-1))/24.0
      enddo
      do j=1,n
         yy(j) = imag(ya(j))
      enddo
      call spline(xa,yy,n,yp1,ypn,y2)
      zi=0.
      do j=2,n
         zi=zi+0.5*(xa(j)-xa(j-1))*(yy(j)+yy(j-1))-
     >      (xa(j)-xa(j-1))**3*(y2(j)+y2(j-1))/24.0
      enddo
      zspline1d = zr + zi*(0.,1.)
      return
      end
