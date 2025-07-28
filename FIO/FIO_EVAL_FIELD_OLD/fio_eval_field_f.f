C     ========================================================
C     read field data ad mapping data
C     YQ Liu, 03/24/2017
C     =======================================================
      subroutine fio_init_field_f(ierr)
      use fio_eval_field_fm

      implicit none
      integer::ierr
      integer::ir,iz,i,j,k
      integer::nsp,nsv
      real*8,dimension(8)::rtmp

      ierr = 0

C     read mapping data from (R,Z) to (S,CHI)
      open(99,file='SCHIMESH_RECTRZ.IN',status='old',form='formatted')
      read(99,*) nr,nz,raxis,zaxis

      allocate( r(nr,nz),z(nr,nz),s(nr,nz),chi(nr,nz) )

      do ir=1,nr
      do iz=1,nz
         read(99,*) r(ir,iz),z(ir,iz),s(ir,iz),chi(ir,iz) 
      enddo
      enddo

      close(99)

C     read B-field data
      open(99,file='BPLASMA_MARSF.IN',status='old',form='formatted')
      read(99,*) n,mmaxe,m1,m2,nsp,nsv
      mmaxp = m2-m1 + 1
      ns    = nsp + nsv
      
      allocate(cs(ns),csm(ns))
      allocate(rmi(ns,mmaxe),zmi(ns,mmaxe),rmm(ns,mmaxe),zmm(ns,mmaxe))
      allocate(b1(ns,mmaxp),b2(ns,mmaxp))

      do i=1,ns
         read(99,*) cs(i),csm(i),rtmp(1)
      enddo

      do j=1,mmaxe
      do i=1,ns
         read(99,*) (rtmp(k),k=1,8)
         rmi(i,j) = cmplx(rtmp(1),rtmp(2))
         zmi(i,j) = cmplx(rtmp(3),rtmp(4))
         rmm(i,j) = cmplx(rtmp(5),rtmp(6))
         zmm(i,j) = cmplx(rtmp(7),rtmp(8))
      enddo
      enddo

      do j=1,mmaxp
      do i=1,ns
         read(99,*) (rtmp(k),k=1,4)
         b1(i,j) = cmplx(rtmp(1),rtmp(2))
         b2(i,j) = cmplx(rtmp(3),rtmp(4))
C        b3(i,j) = cmplx(rtmp(5),rtmp(6))
      enddo
      enddo

      close(99)

      return 
      end

C     ========================================================
C     free memory
C     YQ Liu, 03/24/2017
C     =======================================================
      subroutine fio_finish_field_f
      use fio_eval_field_fm

      deallocate(r,z,s,chi)
      deallocate(cs,csm)
      deallocate(rmi,zmi,rmm,zmm)
      deallocate(b1,b2)

      return
      end

C     ========================================================
C     compute field value at a single point in (R,PHI,Z) space
C     based on MARS-Q data
C     YQ Liu, 03/23/2017
C     =======================================================
      subroutine fio_eval_field_f(type_field,xrphz,brphz,ierr)
      use fio_eval_field_fm

      implicit none
      integer::type_field,ierr,k
      real*8 ::dx,phi0
      real*8,dimension(*)::xrphz,brphz
      real*8,dimension(2,5)::x
      complex*16::ctmp,br,bz,bphi
      complex*16,dimension(2,5)::y

      phi0 = xrphz(2)

      if (type_field.eq.1) then
         do k=1,5
            x(1,k) = xrphz(1)
            x(2,k) = xrphz(3)
         enddo

         dx = 1.e-4

         x(1,2) = x(1,1) - dx
         x(1,3) = x(1,1) + dx

         x(2,4) = x(2,1) - dx
         x(2,5) = x(2,1) + dx
         
         do k=1,5
            call fio_eval_a_f(x(:,k),y(:,k),ierr)
         enddo

         br   = (0.,1.)*n*y(2,1)/x(1,1)
         bz   =-(0.,1.)*n*y(1,1)/x(1,1)
         bphi = (y(1,5)-y(1,4)-y(2,3)+y(2,2))/2./dx

         ctmp = exp(-(0.,1.)*n*phi0)
         brphz(1) = 2.*real(br*ctmp)
         brphz(2) =-2.*real(bphi*ctmp)
         brphz(3) = 2.*real(bz*ctmp)
      endif

      if (type_field.eq.2) then
         x(1,1) = xrphz(1)
         x(2,1) = xrphz(3)

         call fio_eval_a_f(x(:,1),y(:,1),ierr)
         
         ctmp = exp(-(0.,1.)*n*phi0)
         brphz(1) = 2.*real(y(1,1)*ctmp)
         brphz(2) = 0.
         brphz(3) = 2.*real(y(2,1)*ctmp)
      endif

      return
      end

C     ========================================================
C     compute A-field at a single point in (R,Z) space
C     based on MARS-Q data
C     YQ Liu, 03/23/2017
C     =======================================================
      subroutine fio_eval_a_f(xrphz,arphz,ierr)
      use fio_eval_field_fm

      implicit none
      integer::ierr
      integer::ir,iz,i,j,k
      real*8 ::r0,r1,r2,z0,z1,z2,s0,chi0
      real*8 ::drds,dzds,drdc,dzdc,jac,h1,h2,h3,h4,h5,d
      real*8 ::chi00,chi01,chi10,chi11
      real*8,dimension(*)::xrphz
      complex*16,dimension(*)::arphz
      complex*16::y1,y2,y3,y0,yp,ctmp,b10,b20

      ierr = 0

C     for given (r0,z0) find corresponding (s0,chi0)
      r0   = xrphz(1)
      z0   = xrphz(2)

      k  = 0
      do ir=1,nr-1
         if ( (r(ir,1)-r0)*(r(ir+1,1)-r0).le.0. ) k=ir
      enddo
      
      if (k.eq.0) then
         ierr = 1
         write(*,*) 'fio_eval_field_f: ierr=', ierr
         stop
      endif

      r1 = r(k,1)
      r2 = r(k+1,1)
      ir = k

      k  = 0
      do iz=1,nz-1
         if ( (z(1,iz)-z0)*(z(1,iz+1)-z0).le.0. ) k=iz
      enddo
      
      if (k.eq.0) then
         ierr = 2
         write(*,*) 'fio_eval_field_f: ierr=', ierr
         stop
      endif

      z1 = z(1,k)
      z2 = z(1,k+1)
      iz = k

      s0 = ( s(ir,iz)*    (r2-r0)*(z2-z0) +
     &       s(ir,iz+1)*  (r2-r0)*(z0-z1) +
     &       s(ir+1,iz)*  (r0-r1)*(z2-z0) +
     &       s(ir+1,iz+1)*(r0-r1)*(z0-z1) )/(r2-r1)/(z2-z1)

      chi00 = chi(ir,iz)
      chi01 = chi(ir,iz+1)
      chi10 = chi(ir+1,iz)
      chi11 = chi(ir+1,iz+1)

      if ( (max(chi00,chi01,chi10,chi11) -
     &      min(chi00,chi01,chi10,chi11)).lt.3.14 ) then 
     &    
         chi0 = ( chi00*(r2-r0)*(z2-z0) +
     &            chi01*(r2-r0)*(z0-z1) +
     &            chi10*(r0-r1)*(z2-z0) +
     &            chi11*(r0-r1)*(z0-z1) )/(r2-r1)/(z2-z1)
      else
         ctmp = ( exp((0.,1.)*chi00)*(r2-r0)*(z2-z0) +
     &            exp((0.,1.)*chi01)*(r2-r0)*(z0-z1) +
     &            exp((0.,1.)*chi10)*(r0-r1)*(z2-z0) +
     &            exp((0.,1.)*chi11)*(r0-r1)*(z0-z1) )/(r2-r1)/(z2-z1)
         chi0 = datan2(imag(ctmp),real(ctmp))
         if (chi0.lt.0.) chi0 = chi0 + 2.*acos(-1.)
      endif

      if (s0.lt.cs(2)) then
         chi0 = datan2(z0-zaxis,r0-raxis)
         if (chi0.lt.0.) chi0 = chi0 + 2.*acos(-1.)
      endif

C     on the s-mesh, compute dR/ds,dZ/ds,dR/dchi,dZ/dchi, at (s0,chi0)
C     using 3 points: 2 integer  points and 1 middle point
      k = 0
      do i=1,ns-1
         if ( (cs(i)-s0)*(cs(i+1)-s0).le.0. ) k=i
      enddo
      if (k.eq.0) then
         ierr = 3
         write(*,*) 'fio_eval_field_f: ierr=', ierr
         write(*,*) 's0=',s0
         stop
      endif

      i    = k
      h1   = cs(i)   - s0
      h2   = csm(i)  - s0
      h3   = cs(i+1) - s0
      d    = (h1-h2)*(h2-h3)*(h3-h1)
      drds = 0.
      dzds = 0.
      drdc = 0.
      dzdc = 0.

      do j=1,mmaxe
         y1=rmi(i,j)
         y2=rmm(i,j)
         y3=rmi(i+1,j)
         y0=(h3*h2*(h3-h2)*y1+h1*h3*(h1-h3)*y2+h2*h1*(h2-h1)*y3)/d
         yp=((h2-h3)*(h2+h3)*y1+(h3-h1)*(h3+h1)*y2+(h1-h2)*(h1+h2)*y3)/d
         
         if (j.eq.1) then
            drds = drds + real(yp)
         else
            ctmp = exp((0.,1.)*(j-1)*chi0)
            drds = drds + 2.*real(yp*ctmp)
            drdc = drdc + 2.*real(y0*(0.,1.)*(j-1)*ctmp)
         endif

         y1=zmi(i,j)
         y2=zmm(i,j)
         y3=zmi(i+1,j)
         y0=(h3*h2*(h3-h2)*y1+h1*h3*(h1-h3)*y2+h2*h1*(h2-h1)*y3)/d
         yp=((h2-h3)*(h2+h3)*y1+(h3-h1)*(h3+h1)*y2+(h1-h2)*(h1+h2)*y3)/d
         
         if (j.eq.1) then
            dzds = dzds + real(yp)
         else
            dzds = dzds + 2.*real(yp*ctmp)
            dzdc = dzdc + 2.*real(y0*(0.,1.)*(j-1)*ctmp)
         endif
      enddo

C     compute jacobian
      jac = drds*dzdc - drdc*dzds

C     compute b1,b2,b3 at (s0,chi0)
      if (s0.lt.1..and.i.gt.1) i=i-1 
      h1   = cs(i)    - s0
      h2   = cs(i+1)  - s0
      h3   = cs(i+2)  - s0
      h4   = csm(i)   - s0
      h5   = csm(i+1) - s0
      d    = (h1-h2)*(h2-h3)*(h3-h1)

      b10 = (0.,0.)
      b20 = (0.,0.)      
      do j=1,mmaxp
         ctmp = exp((0.,1.)*(j-1+m1)*chi0)

         y1 = b1(i,j)
         y2 = b1(i+1,j)
         y3 = b1(i+2,j)
         y0=(h3*h2*(h3-h2)*y1+h1*h3*(h1-h3)*y2+h2*h1*(h2-h1)*y3)/d
C        yp=((h2-h3)*(h2+h3)*u1+(h3-h1)*(h3+h1)*u2+(h1-h2)*(h1+h2)*u3)/d
         b10 = b10 + y0*ctmp

         y1 = (h5*b2(i,j)-h4*b2(i+1,j))/(h5-h4)
         b20 = b20 + y1*ctmp
      enddo
      
C     compute Ar,Az at xrphiz=(r0,z0)
      arphz(1) = (b10*dzds+b20*dzdc)/jac/(-(0.,1.)*n)
      arphz(2) = (b10*drds+b20*drdc)/jac/((0.,1.)*n)   
      
      return
      end

