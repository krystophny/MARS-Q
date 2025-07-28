	Module NTV_lib
C******************************************************************
!	Some subroutines needed for NTV calculation
!
!           some of them comes from the book 
!              Numerical recipes in Fortran
!              Second Edition
!           Press Syndicate of the University of Cambridge
!       
!       by Y. Sun
C******************************************************************
	implicit none

	contains

	function cumsum(n,arr)
C***************************************************
C        cumulative sum on an array
C***************************************************
	integer n,j
	real arr(n),cumsum(n)
        cumsum(1)=arr(1)
	do j=2,n
	   cumsum(j)=cumsum(j-1)+arr(j)
	end do
	end function cumsum

	subroutine trans_RZ(fc,fs,nr,M0,n_chi,chi,fr)
C***************************************************
C     reconstruct the (R,Z)profile from MAS-F output
C             to real space
C          2D (r,chi)
C***************************************************
	integer M0,nr,n_chi,i,j
	real fc(nr,M0),fs(nr,M0),chi(n_chi),fr(nr,n_chi)
        
	do 10 i=1,nr
	     fr(i,:)=0.0     
	   do 20 j=2,M0
	     fr(i,:)=fr(i,:)+2.0*(fc(i,j)*cos((j-1.0)*chi)-fs(i,j)*
     &               sin((j-1.0)*chi))
 20	   continue
	     fr(i,:)=fc(i,1)+fr(i,:)
 10	continue
	end subroutine trans_RZ

	subroutine trans_2d(fc,fs,M0,M,n,n_chi,n_phi,chi,phi,fr,fi)
C***************************************************
C     reconstruct the Xi profile from MAS-F output
C             to real space
C          from f(M)
C          to f(chi,phi) on each flux surface
C***************************************************
	integer M0,n,nr,n_chi,n_phi,i,j,M(M0),k
	real fc(M0),fs(M0),chi(n_chi),phi(n_chi,n_phi)
	real fr(n_chi,n_phi),fi(n_chi,n_phi)
c$$$          print*,'fc',fc
c$$$	  print*,'M',M
c$$$         print*,'chi',chi
c$$$	  print*,'n',n
c$$$	  print*,'phi',phi(:,1)
      do 10 i=1,n_phi  
	 do 10 k=1,n_chi
            fr(k,i)=0.0
            fi(k,i)=0.0
	   do 10 j=1,M0
	    fr(k,i)=fr(k,i)+fc(j)*cos(M(j)*chi(k)+n*phi(k,i))-
     .                      fs(j)*sin(M(j)*chi(k)+n*phi(k,i))
	    fi(k,i)=fi(k,i)+fs(j)*cos(M(j)*chi(k)+n*phi(k,i))+
     .                      fc(j)*sin(M(j)*chi(k)+n*phi(k,i))
       !    print*,j,M(j),fc(j),fs(j),fr(k,i)
 10   continue
        !   print*,'fr',fr(:,1)
        !   pause
	end subroutine trans_2d

	subroutine trans_1d(fc,fs,M0,M,n,n_chi,chi,dphi,fr,fi)
C***************************************************
C     reconstruct the Xi profile from MAS-F output
C             to real space
C          from f(M)
C          to f(chi) on each flux surface
C***************************************************
	integer M0,n,n_chi,j,M(M0),k
	real fc(M0),fs(M0),chi(n_chi),dphi(n_chi)
	real fr(n_chi),fi(n_chi)
	do 10 k=1,n_chi
            fr(k)=0.0
            fi(k)=0.0
	   do 10 j=1,M0
	    fr(k)=fr(k)+fc(j)*cos(M(j)*chi(k)-n*dphi(k))-
     .                      fs(j)*sin(M(j)*chi(k)-n*dphi(k))
	    fi(k)=fi(k)+fs(j)*cos(M(j)*chi(k)-n*dphi(k))+
     .                      fc(j)*sin(M(j)*chi(k)-n*dphi(k))
 10     continue
	end subroutine trans_1d

      subroutine calint(f,x,nx,intf)
C***************************************************
C     calculate \int (fdx)  for non-uniform x
C***************************************************
      implicit none
      integer nx,i
      real x(nx),dx(nx-1)
      complex*16 f(nx),intf,f1(nx-1)
      f1=0.5*(f(1:(nx-1))+f(2:nx))
      dx=abs(x(2:nx)-x(1:(nx-1)))
      intf=sum(f1*dx)
      end subroutine calint

      subroutine cald(f,nx,ny,x,dfdx)
C***************************************************
C     calculate df(x,y)/dx  for non-uniform x
C***************************************************
      implicit none
      integer nx,ny,i
      real x(nx),f(nx,ny),dxu(nx),dxd(nx)
      real dfu(nx,ny),dfd(nx,ny),dfdx(nx,ny)
      do i=2,nx
         dxd(i)=x(i)-x(i-1)
         dxu(i-1)=dxd(i)
         dfd(i,:)=f(i,:)-f(i-1,:)
         dfu(i-1,:)=dfd(i,:)
      end do
      dxu(nx)=dxu(nx-1)
      dxd(1)=dxd(2)
      dfu(nx,:)=dfu(nx-1,:)
      dfd(1,:)=dfd(2,:)
      do i=1,ny
         dfdx(:,i)=(dxd*dfu(:,i)/dxu+dxu*dfd(:,i)/dxd)/(dxu+dxd)
      end do
      end subroutine cald

      subroutine cald_2d(f,nx,ny,x,y,dfdx,dfdy)
C***************************************************
C     calculate df(x,y)/dx  for non-uniform x
C               df(x,y)/dy  for non-uniform y
C***************************************************
      implicit none
      integer nx,ny,i
      real x(nx),y(ny),f(nx,ny),dxu(nx),dxd(nx)
      real dfu(nx,ny),dfd(nx,ny),dfdx(nx,ny)
      real dyl(ny),dyr(ny)
      real dfl(nx,ny),dfr(nx,ny),dfdy(nx,ny)
      do i=2,nx
         dxd(i)=x(i)-x(i-1)
         dxu(i-1)=dxd(i)
         dfd(i,:)=f(i,:)-f(i-1,:)
         dfu(i-1,:)=dfd(i,:)
      end do
      dxu(nx)=dxu(nx-1)
      dxd(1)=dxd(2)
      dfu(nx,:)=dfu(nx-1,:)
      dfd(1,:)=dfd(2,:)
      do i=1,ny
         dfdx(:,i)=(dxd*dfu(:,i)/dxu+dxu*dfd(:,i)/dxd)/(dxu+dxd)
      end do

      do i=2,ny
         dyl(i)=y(i)-y(i-1)
         dyr(i-1)=dyl(i)
         dfl(:,i)=f(:,i)-f(:,i-1)
         dfr(:,i-1)=dfl(:,i)
      end do
      dyr(ny)=dyr(ny-1)
      dyl(1)=dyl(2)
      dfr(:,ny)=dfr(:,ny-1)
      dfl(:,1)=dfl(:,2)
      do i=1,nx
         dfdy(i,:)=(dyl*dfr(i,:)/dyr+dyr*dfl(i,:)/dyl)/(dyr+dyl)
      end do

      end subroutine cald_2d

      subroutine cald_1d(f,nx,x,dfdx)
C***************************************************
C     calculate df(x)/dx  for non-uniform x
C***************************************************
      implicit none
      integer nx,i
      real x(nx),f(nx),dxu(nx),dxd(nx)
      real dfu(nx),dfd(nx),dfdx(nx)
      do i=2,nx
         dxd(i)=x(i)-x(i-1)
         dxu(i-1)=dxd(i)
         dfd(i)=f(i)-f(i-1)
         dfu(i-1)=dfd(i)
      end do
      dxu(nx)=dxu(nx-1)
      dxd(1)=dxd(2)
      dfu(nx)=dfu(nx-1)
      dfd(1)=dfd(2)
         dfdx=(dxd*dfu/dxu+dxu*dfd/dxd)/(dxu+dxd)
      end subroutine cald_1d

      subroutine xgrid(nx,x_min,x_max,x)
C***************************************************
C      generate equally spaced x grid
C        x_min < x(1:n) < x_max
C***************************************************
      implicit none
      integer i,nx
      real x_max,x_min,x(nx),dx
      dx=(x_max-x_min)/nx
      do i=1,nx
         x(i)=(i-0.5)*dx+x_min
      end do
      end subroutine xgrid

      subroutine xgrid_f(nx,x_min,x_max,x)
C***************************************************
C      generate equally spaced x grid
C        x_min <= x(1:n) <= x_max
C***************************************************
      implicit none
      integer i,nx
      real x_max,x_min,x(nx),dx
      dx=(x_max-x_min)/(nx-1.0)
      do i=1,nx
         x(i)=(i-1.0)*dx+x_min
      end do
      end subroutine xgrid_f

	subroutine fft1(n,fx,fk)
C**********************************************
C       1-D FFT  
C       use four1      
C            by Sun          
C**********************************************
	integer n,i,j
        real fx(n),ndata(2*n)
	complex fk(n)
	ndata=0.0
	do i=1,n
	   j=2*i-1
	   ndata(j)=fx(i)
	end do
	call four1(ndata,n,1)
	do i=1,n
	   j=2*i-1
	   fk(i)=cmplx(ndata(j),ndata(j+1))
	end do
	end subroutine fft1

	subroutine fftn(ndim,n,fx,fk)
C**********************************************
C     ndim-D FFT  
C    use fourn   
C            by Sun          
C**********************************************
	complex fk(*)
        real fx(*)
	integer ndim,n(ndim)
	integer i,j,ntot
	real,allocatable::ndata(:)
	 ntot=1
	 do i=1,ndim
          ntot=ntot*n(i)
	 end do
	 ntot=2*ntot
    !     print*, n,ntot
	 allocate(ndata(ntot))

          j=0
	do 10 i=1,ntot/2
             j=2*i-1
	    ndata(j)=fx(i)
	    ndata(j+1)=0.0
 10	continue
	call fourn(ndata,n,ndim,1)
           j=0
	do 20 i=1,ntot/2
                j=2*i-1
             fk(i)=cmplx(ndata(j),ndata(j+1))
 20	continue
	deallocate(ndata)
	end subroutine fftn

	FUNCTION erf(x)
	REAL erf,x
C******************************************************************
!	C USES gammp
!	Returns the error function erf(x).
C******************************************************************
!	REAL gammp	
	if(x.lt.0.)then
	erf=-gammp(.5,x**2)
	else
	erf=gammp(.5,x**2)
	endif
	return
	END FUNCTION erf

	subroutine ellipke(k2,nk,K,E)
C******************************************************************
!	calculate complete elliptic intergral of 
!      the first and second kinds ( K and E )
!        note: the argument k2=kappa**2
!            k2=k2(1:nk)
!     use ellf elle
!       by Y. Sun
C******************************************************************
	real k2(*),K(*),E(*)
	integer i,nk 
	real k0
	do i=1,nk
	   k0=sqrt(k2(i))	!  change the argument to k^2
	   K(i)=ellf(3.1415926535897932/2.0,k0)
	   E(i)=elle(3.1415926535897932/2.0,k0)
	end do
	end subroutine ellipke

	FUNCTION ellf(phi,ak)
	REAL ellf,ak,phi
C******************************************************************
!	C USES rf
!	Legendre elliptic integral of the 1st kind F(.; k), evaluated using Carlson's function RF .
!	The argument ranges are 0 <=phi<=pi/2, 0 <=k sin(phi)<=1.
C******************************************************************
!	REAL s,rf     ! for internal function rf, it is not necessary to define the data type
	REAL s
	s=sin(phi)
	ellf=s*rf(cos(phi)**2,(1.-s*ak)*(1.+s*ak),1.)
	return
	END FUNCTION ellf

	FUNCTION elle(phi,ak)
	REAL elle,ak,phi
C******************************************************************
!	C USES rd,rf
!	Legendre elliptic integral of the 2nd kind E(.; k), evaluated using Carlson's functions RD
!	and RF . The argument ranges are  0 <=phi<=pi/2, 0 <=k sin(phi)<=1.
C******************************************************************
!	REAL cc,q,s,rd,rf
      REAL cc,q,s

	s=sin(phi)
	cc=cos(phi)**2
	q=(1.-s*ak)*(1.+s*ak)
	elle=s*(rf(cc,q,1.)-((s*ak)**2)*rd(cc,q,1.)/3.)
	return
	END FUNCTION elle

	FUNCTION rd(x,y,z)
C******************************************************************
!	Computes Carlson's elliptic integral of the second kind, RD(x; y; z). x and y must be
!	nonnegative, and at most one can be zero. z must be positive. TINY must be at least twice
!	the negative 2/3 power of the machine overflow limit. BIG must be at most 0:1.ERRTOL
!	times the negative 2/3 power of the machine underflow limit.
C******************************************************************
	REAL rd,x,y,z,ERRTOL,TINY,BIG,C1,C2,C3,C4,C5,C6
	PARAMETER (ERRTOL=.05,TINY=1.e-25,BIG=4.5E21,C1=3./14.,C2=1./6.,
     .   C3=9./22.,C4=3./26.,C5=.25*C3,C6=1.5*C4)

	REAL alamb,ave,delx,dely,delz,ea,eb,ec,ed,ee,fac,sqrtx,sqrty,
     .     sqrtz,sum,xt,yt,zt
	if(min(x,y).lt.0..or.min(x+y,z).lt.TINY.or.
     .    max(x,y,z).gt.BIG) stop 'invalid arguments in rd'
	xt=x
	yt=y
	zt=z
	sum=0.
	fac=1.
1     continue
	sqrtx=sqrt(xt)
	sqrty=sqrt(yt)
	sqrtz=sqrt(zt)
	alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
	sum=sum+fac/(sqrtz*(zt+alamb))
	fac=.25*fac
	xt=.25*(xt+alamb)
	yt=.25*(yt+alamb)
	zt=.25*(zt+alamb)
	ave=.2*(xt+yt+3.*zt)
	delx=(ave-xt)/ave
	dely=(ave-yt)/ave
	delz=(ave-zt)/ave
	if(max(abs(delx),abs(dely),abs(delz)).gt.ERRTOL)goto 1
	ea=delx*dely
	eb=delz*delz
	ec=ea-eb
	ed=ea-6.*eb
	ee=ed+ec+ec
	rd=3.*sum+fac*(1.+ed*(-C1+C5*ed-C6*delz*ee)
     . +delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*sqrt(ave))
	return
	END FUNCTION rd

	FUNCTION rf(x,y,z)
C******************************************************************
!	Computes Carlson's elliptic integral of the rst kind, RF (x; y; z). x, y, and z must be
!	nonnegative, and at most one can be zero. TINY must be at least 5 times the machine
!	underflow limit, BIG at most one fth the machine overflow limit.
C******************************************************************
	REAL rf,x,y,z,ERRTOL,TINY,BIG,THIRD,C1,C2,C3,C4
	PARAMETER (ERRTOL=.08,TINY=1.5e-38,BIG=3.E37,THIRD=1./3.,
     .	C1=1./24.,C2=.1,C3=3./44.,C4=1./14.)
	REAL alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt
	if(min(x,y,z).lt.0..or.min(x+y,x+z,y+z).lt.TINY.or.
     . max(x,y,z).gt.BIG) stop 'invalid arguments in rf'
	xt=x
	yt=y
	zt=z
1     continue
	sqrtx=sqrt(xt)
	sqrty=sqrt(yt)
	sqrtz=sqrt(zt)
	alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
	xt=.25*(xt+alamb)
	yt=.25*(yt+alamb)
	zt=.25*(zt+alamb)
	ave=THIRD*(xt+yt+zt)
	delx=(ave-xt)/ave
	dely=(ave-yt)/ave
	delz=(ave-zt)/ave
	if(max(abs(delx),abs(dely),abs(delz)).gt.ERRTOL)goto 1
	e2=delx*dely-delz**2
	e3=delx*dely*delz
	rf=(1.+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave)
	return
	END FUNCTION rf

	FUNCTION gammp(a,x)
C******************************************************************
!C USES gcf,gser
!Returns the incomplete gamma function P(a; x).
C******************************************************************
	REAL a,gammp,x
	REAL gammcf,gamser,gln
	if(x.lt.0..or.a.le.0.) stop 'bad arguments in gammp'
	if(x.lt.a+1.)then !Use the series representation.
	call gser(gamser,a,x,gln)
	gammp=gamser
	else !Use the continued fraction representation
	call gcf(gammcf,a,x,gln)
	gammp=1.-gammcf !and take its complement.
	endif
	return
	END FUNCTION gammp

	SUBROUTINE gser(gamser,a,x,gln)
	INTEGER ITMAX
	REAL a,gamser,gln,x,EPS
	PARAMETER (ITMAX=100,EPS=3.e-7)
C******************************************************************
!	C USES gammln
!	Returns the incomplete gamma function P(a; x) evaluated by its series representation as
!	gamser. Also returns ln..(a) as gln.
C******************************************************************
	INTEGER n
!	REAL ap,del,sum,gammln
      REAL ap,del,sum
	gln=gammln(a)
	if(x.le.0.)then
	if(x.lt.0.) stop 'x < 0 in gser'
	gamser=0.
	return
	endif
	ap=a
	sum=1./a
	del=sum
	do 11 n=1,ITMAX
	ap=ap+1.
	del=del*x/ap
	sum=sum+del
	if(abs(del).lt.abs(sum)*EPS)goto 1
11    continue
	stop 'a too large, ITMAX too small in gser'
1     gamser=sum*exp(-x+a*log(x)-gln)
	return
	END SUBROUTINE gser


	SUBROUTINE gcf(gammcf,a,x,gln)
	INTEGER ITMAX
	REAL a,gammcf,gln,x,EPS,FPMIN
	PARAMETER (ITMAX=100,EPS=3.e-7,FPMIN=1.e-30)
C******************************************************************
!	C USES gammln
!	Returns the incomplete gamma function Q(a; x) evaluated by its continued fraction representation
!	as gammcf. Also returns ln..(a) as gln.
!	Parameters: ITMAX is the maximum allowed number of iterations; EPS is the relative accuracy;
!	FPMIN is a number near the smallest representable floating-point number.
C******************************************************************
	INTEGER i
!	REAL an,b,c,d,del,h,gammln
      REAL an,b,c,d,del,h
	gln=gammln(a)
	b=x+1.-a !Set up for evaluating continued fraction by modied
	c=1./FPMIN !Lentz's method (x5.2) with b0 = 0.
	d=1./b
	h=d
	do 11 i=1,ITMAX !Iterate to convergence.
	an=-i*(i-a)
	b=b+2.
	d=an*d+b
	if(abs(d).lt.FPMIN)d=FPMIN
	c=b+an/c
	if(abs(c).lt.FPMIN)c=FPMIN
	d=1./d
	del=d*c
	h=h*del
	if(abs(del-1.).lt.EPS)goto 1
11    continue
	stop 'a too large, ITMAX too small in gcf'
1     gammcf=exp(-x+a*log(x)-gln)*h !Put factors in front.
	return
	END SUBROUTINE gcf


	FUNCTION gammln(xx)
	REAL gammln,xx
C******************************************************************
!	Returns the value ln[Gamma(xx)] for xx > 0.
C******************************************************************
	INTEGER j
	DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
!	Internal arithmetic will be done in double precision, a nicety that you can omit if 
!five-figure accuracy is good enough.
	SAVE cof,stp
	DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     * 24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     * -.5395239384953d-5,2.5066282746310005d0/
	x=xx
	y=x
	tmp=x+5.5d0
	tmp=(x+0.5d0)*log(tmp)-tmp
	ser=1.000000000190015d0
	do 11 j=1,6
	y=y+1.d0
	ser=ser+cof(j)/y
11    continue
	gammln=tmp+log(stp*ser/x)
	return
	END FUNCTION gammln

	SUBROUTINE fourn(data,nn,ndim,isign)
	INTEGER isign,ndim,nn(ndim)
	REAL data(*)
C       Replaces data by its ndim-dimensional discrete Fourier transform, if isign is input as
C       1. nn(1:ndim) is an integer array containing the lengths of each dimension (number of
C       complex values), which MUST all be powers of 2. data is a real array of length twice the
C       product of these lengths, in which the data are stored as in a multidimensional complex
C       FORTRAN array. If isign is input as −1, data is replaced by its inverse transform times
C       the product of the lengths of all dimensions.
	INTEGER i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1,ifp2,ip1,ip2,
     .    ip3,k1,k2,n,nprev,nrem,ntot
	REAL tempi,tempr
!	COMPLEX  tempi,tempr
	DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp    !Double precision for trigonometric recurrences.
	ntot=1
                           
	do 11 idim=1,ndim !Compute total number of complex values.
	   ntot=ntot*nn(idim)
 11	continue
	nprev=1
	do 18 idim=1,ndim !Main loop over the dimensions.
	   n=nn(idim)
	   nrem=ntot/(n*nprev)
	   ip1=2*nprev
	   ip2=ip1*n
	   ip3=ip2*nrem
	   i2rev=1
	   do 14 i2=1,ip2,ip1 !This is the bit-reversal section of the routine.
	      if(i2.lt.i2rev)then
		 do 13 i1=i2,i2+ip1-2,2
		    do 12 i3=i1,ip3,ip2
		       i3rev=i2rev+i3-i2
		       tempr=data(i3)
		       tempi=data(i3+1)
		       data(i3)=data(i3rev)
		       data(i3+1)=data(i3rev+1)
		       data(i3rev)=tempr
		       data(i3rev+1)=tempi
 12		    continue
 13		 continue
	      endif
	      ibit=ip2/2
 1	      if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
		 i2rev=i2rev-ibit
		 ibit=ibit/2
		 goto 1
	      endif
	      i2rev=i2rev+ibit
 14	   continue
	   ifp1=ip1 !Here begins the Danielson-Lanczos section of the routine.
 2	   if(ifp1.lt.ip2)then
	      ifp2=2*ifp1
	      theta=isign*6.28318530717959d0/(ifp2/ip1) !Initialize for the trig. recurrence.
	      wpr=-2.d0*sin(0.5d0*theta)**2 
	      wpi=sin(theta)
	      wr=1.d0
	      wi=0.d0
	      do 17 i3=1,ifp1,ip1
		 do 16 i1=i3,i3+ip1-2,2
		    do 15 i2=i1,ip3,ifp2
		       k1=i2               !Danielson-Lanczos formula:
		       k2=k1+ifp1
		       tempr=sngl(wr)*data(k2)-sngl(wi)*data(k2+1)
		       tempi=sngl(wr)*data(k2+1)+sngl(wi)*data(k2)
		       data(k2)=data(k1)-tempr
		       data(k2+1)=data(k1+1)-tempi

		       data(k1)=data(k1)+tempr
		       data(k1+1)=data(k1+1)+tempi
 15		    continue
 16		 continue
		 wtemp=wr !Trigonometric recurrence.
		 wr=wr*wpr-wi*wpi+wr
		 wi=wi*wpr+wtemp*wpi+wi
 17	      continue
	      ifp1=ifp2
	      goto 2
	   endif
	   nprev=n*nprev
 18	continue
	return
	END SUBROUTINE fourn

	SUBROUTINE four1(data,nn,isign)
	INTEGER isign,nn
	REAL data(2*nn)
C	Replaces data(1:2*nn) by its discrete Fourier transform, if isign is input as 1; or replaces
C       data(1:2*nn) by nn times its inverse discrete Fourier transform, if isign is input as −1.
C       data is a complex array of length nn or, equivalently, a real array of length 2*nn. nn
C	MUST be an integer power of 2 (this is not checked for !).
	INTEGER i,istep,j,m,mmax,n
	REAL tempi,tempr
        DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp 
                  ! Double precision for the trigonometric recurrences.
	n= 2*nn
	j=1
	do 11 i=1,n,2		!This is the bit-reversal section of the routine.
	   if(j.gt.i)then
	      tempr=data(j)	!Exchange the two complex numbers.
	      tempi=data(j+1)
	      data(j)=data(i)
	      data(j+1)=data(i+1)
	      data(i)=tempr
	      data(i+1)=tempi
	   endif
	   m=n/2
 1	   if ((m.ge.2).and.(j.gt.m)) then
	      j=j-m
	      m=m/2
	      goto 1
	   endif
	   j=j+m
 11	continue

	mmax=2			!Here begins the Danielson-Lanczos section of the routine.
 2	if (n.gt.mmax) then	!Outer loop executed log2 nn times.
	istep=2*mmax
	theta=6.28318530717959d0/(isign*mmax) !Initialize for the trigonometric recurrence.
	wpr=-2.d0*sin(0.5d0*theta)**2 
	wpi=sin(theta)
	wr=1.d0
	wi=0.d0
	do 13 m=1,mmax,2	!Here are the two nested inner loops.
	   do 12 i=m,n,istep
	      j=i+mmax		!This is the Danielson-Lanczos formula:
	      tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
	      tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
	      data(j)=data(i)-tempr
	      data(j+1)=data(i+1)-tempi
	      data(i)=data(i)+tempr
	      data(i+1)=data(i+1)+tempi
 12	continue
	wtemp=wr		!Trigonometric recurrence.
	wr=wr*wpr-wi*wpi+wr
	wi=wi*wpr+wtemp*wpi+wi
 13	continue
	mmax=istep
	goto 2			!Not yet done.
	endif			!All done.
	return
	END SUBROUTINE four1


	end



