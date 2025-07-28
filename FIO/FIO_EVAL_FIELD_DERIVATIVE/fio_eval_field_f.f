C     ========================================================
C     read field data and mapping data
C     YQ Liu, 03/24/2017
C     =======================================================
      subroutine fio_init_field_f(ierr)
      use fio_eval_field_fm

      implicit none
      integer::ierr
      integer::ir,iz,i,j,k,n
      integer::nsp,nsv,mmaxpm
      real*8,dimension(8)::rtmp
      character(len=128) filename

      ierr     = 0
      pi_value = acos(-1.)
      nmax     = 1

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

      allocate(nn(nmax),ns(nmax),mmaxp(nmax),m1(nmax),m2(nmax))

C     first read B-field data to get dimensions 
      do n=1,nmax
         write(filename,"(A16,I0.1,A3)") "BPLASMA_MARSF_",n,".IN"
         open(99,file=trim(filename),status='old',form='formatted')
         read(99,*) nn(n),mmaxe,m1(n),m2(n),nsp,nsv
         close(99)
         mmaxp(n) = m2(n)-m1(n) + 1
         ns(n)    = nsp + nsv
      enddo

C     find max(mmaxp)
      mmaxpm = mmaxp(1)
      nsm    = ns(1)
      do n=1,nmax
         if (mmaxpm.lt.mmaxp(n)) mmaxpm = mmaxp(n)
         if (nsm.lt.ns(n))       nsm = ns(n)
      enddo

C     allocate new arrays
      allocate(cs(nsm,nmax),csm(nsm,nmax))
      allocate(rmi(nsm,mmaxe,nmax),zmi(nsm,mmaxe,nmax),
     &         rmm(nsm,mmaxe,nmax),zmm(nsm,mmaxe,nmax))
      allocate(b1(nsm,mmaxpm,nmax),b2(nsm,mmaxpm,nmax),
     &         b3(nsm,mmaxpm,nmax))
      allocate(x1(nsm,mmaxpm,nmax))

      cs  = 0.
      csm = 0.
      rmi = 0.
      zmi = 0.
      rmm = 0.
      zmm = 0.
      b1  = 0.
      b2  = 0.
      b3  = 0.
      x1  = 0.

C     read B-field data for all toroidal harmonics
C     followed by X1 data for normal displacement of the plasma      
      do n=1,nmax

      write(filename,"(A16,I0.1,A3)") "BPLASMA_MARSF_n",n,".IN"
      open(99,file=trim(filename),status='old',form='formatted')

      read(99,*) nn(n),mmaxe,m1(n),m2(n),nsp,nsv
      do i=1,ns(n)
         read(99,*) cs(i,n),csm(i,n),rtmp(1)
      enddo

      do j=1,mmaxe
      do i=1,ns(n)
         read(99,*) (rtmp(k),k=1,8)
         rmi(i,j,n) = cmplx(rtmp(1),rtmp(2))
         zmi(i,j,n) = cmplx(rtmp(3),rtmp(4))
         rmm(i,j,n) = cmplx(rtmp(5),rtmp(6))
         zmm(i,j,n) = cmplx(rtmp(7),rtmp(8))
      enddo
      enddo

      do j=1,mmaxp(n)
      do i=1,ns(n)
         read(99,*) (rtmp(k),k=1,6)
         b1(i,j,n) = cmplx(rtmp(1),rtmp(2))
         b2(i,j,n) = cmplx(rtmp(3),rtmp(4))
         b3(i,j,n) = cmplx(rtmp(5),rtmp(6))
      enddo
      enddo

      do j=1,mmaxp(n)
      do i=1,nsp+1
         read(99,*) (rtmp(k),k=1,2)
         x1(i,j,n) = cmplx(rtmp(1),rtmp(2))
      enddo
      enddo

      close(99)
      enddo

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
      deallocate(nn,mmaxp,m1,m2)
      deallocate(b1,b2,b3,x1)

      return
      end

C     ========================================================
C     compute the field at a single point in (R,Z) space
C     based on MARS-Q data
C     YQ Liu, 03/23/2017
C     =======================================================
      subroutine fio_eval_field_f(type_field,xrphz,brphz,dbrphz,ierr)
      use fio_eval_field_fm

      implicit none
      integer::type_field,ierr
      complex*16,dimension(3)  ::xrphz,brphz
      complex*16,dimension(3,3)::dbrphz

      integer::ir,iz,i,j,k,n
      real*8 ::r0,r1,r2,z0,z1,z2,s0,chi0,phi0,r0n,z0n
      real*8 ::h1,h2,h3,h4,d0,d1,d2,
     &         h01,h02,h03,h11,h12,h13,h21,h22,h23,h24,
     &         g01,g02,g03,g11,g12,g13
      real*8 ::jac,djacds,djacdc,
     &         drds,dzds,drdc,dzdc,
     &         drds2,drdsc,drdc2,dzds2,dzdsc,dzdc2
      real*8 ::chi00,chi01,chi10,chi11
      complex*16::f1,f2,f3,f4,q0,q1,q2
      complex*16::br,bp,bz,
     &            dbrds,dbrdc,dbrdr,dbrdz,dbrdp,
     &            dbzds,dbzdc,dbzdr,dbzdz,dbzdp,
     &            dbpds,dbpdc,dbpdr,dbpdz,dbpdp
      complex*16::y1,y2,y3,y0,yp,ctmp,b10,b20,b30,
     &            db10ds,db10dc,db10dp,
     &            db20ds,db20dc,db20dp,
     &            db30ds,db30dc,db30dp

      ierr = 0

C     for given (r0,z0) find corresponding (s0,chi0)
      r0   = xrphz(1)
      phi0 = xrphz(2)
      z0   = xrphz(3)
      

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
     &      min(chi00,chi01,chi10,chi11)).lt.pi_value ) then 
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
         if (chi0.lt.0.) chi0 = chi0 + 2.*pi_value
      endif

      if (s0.lt.maxval(cs(2,:))) then
         chi0 = datan2(z0-zaxis,r0-raxis)
         if (chi0.lt.0.) chi0 = chi0 + 2.*pi_value
      endif

      brphz  = 0.
      dbrphz = 0.

C     go through all toroidal harmonics
      do n=1,nmax

C     on the s-mesh, compute the following equilibrium quantities at (s0,chi0)
C     a) jacobian
C     b) first order derivatives of (R,Z,jacobian) w.r.t. (s,chi)
C     c) second order derivatives of (R,Z) w.r.t. (s,chi)

      k = 0
      do i=1,ns(n)-1
         if ( (cs(i,n)-s0)*(cs(i+1,n)-s0).le.0. ) k=i
      enddo
      if (k.eq.0) then
         ierr = 3
         write(*,*) 'fio_eval_field_f: ierr=', ierr
         write(*,*) 's0=',s0
         stop
      endif

      i    = k
      h1   = cs(i,n)    - s0
      h2   = csm(i,n)   - s0
      h3   = cs(i+1,n)  - s0
      h4   = csm(i+1,n) - s0

      d0   = (h1-h2)*(h2-h3)*(h3-h1)
      d1   = d0
      d2   = 0.5*(h4-h3)*(h4-h2)*(h4-h1)*(h3-h2)*(h3-h1)*(h2-h1)
      h01  = h3*h2*(h3-h2)/d0 
      h02  = h1*h3*(h1-h3)/d0
      h03  = h2*h1*(h2-h1)/d0
      h11  = (h2-h3)*(h2+h3)/d1
      h12  = (h3-h1)*(h3+h1)/d1
      h13  = (h1-h2)*(h1+h2)/d1
      h21  = (h2+h3+h4)*(h4-h3)*(h4-h2)*(h3-h2)/d2
      h22  =-(h1+h3+h4)*(h4-h3)*(h4-h1)*(h3-h1)/d2
      h23  = (h1+h2+h4)*(h4-h2)*(h4-h1)*(h2-h1)/d2
      h24  =-(h1+h2+h3)*(h3-h2)*(h3-h1)*(h2-h1)/d2

C     first and second oder derivatives of (R,Z) w.r.t. (s,chi)
      r0n   = 0.
      z0n   = 0.
      drds  = 0.
      dzds  = 0.
      drdc  = 0.
      dzdc  = 0.
      drds2 = 0.
      drdsc = 0.
      drdc2 = 0.
      dzds2 = 0.
      dzdsc = 0.
      dzdc2 = 0.

      do j=1,mmaxe
         f1 = rmi(i,j,n)
         f2 = rmm(i,j,n)
         f3 = rmi(i+1,j,n)
         f4 = rmm(i+1,j,n)
         q0 = h01*f1 + h02*f2 + h03*f3
         q1 = h11*f1 + h12*f2 + h13*f3
         q2 = h21*f1 + h22*f2 + h23*f3 + h24*f4
         
         if (j.gt.1) ctmp = exp((0.,1.)*(j-1)*chi0)

         if (j.eq.1) then
            r0n   = r0n   + real(q0)
            drds  = drds  + real(q1)
            drds2 = drds2 + real(q2)
         else
            r0n   = r0n   + 2.*real(q0*ctmp)
            drds  = drds  + 2.*real(q1*ctmp)
            drdc  = drdc  + 2.*real(q0*(0.,1.)*(j-1)*ctmp)
            drds2 = drds2 + 2.*real(q2*ctmp)
            drdsc = drdsc + 2.*real(q1*(0.,1.)*(j-1)*ctmp)
            drdc2 = drdc2 - 2.*real(q0*(j-1)**2*ctmp)
         endif

         f1 = zmi(i,j,n)
         f2 = zmm(i,j,n)
         f3 = zmi(i+1,j,n)
         f4 = zmm(i+1,j,n)
         q0 = h01*f1 + h02*f2 + h03*f3
         q1 = h11*f1 + h12*f2 + h13*f3
         q2 = h21*f1 + h22*f2 + h23*f3 + h24*f4
         
         if (j.eq.1) then
            z0n   = z0n   + real(q0)
            dzds  = dzds  + real(q1)
            dzds2 = dzds2 + real(q2)
         else
            z0n   = z0n   + 2.*real(q0*ctmp)
            dzds  = dzds  + 2.*real(q1*ctmp)
            dzdc  = dzdc  + 2.*real(q0*(0.,1.)*(j-1)*ctmp)
            dzds2 = dzds2 + 2.*real(q2*ctmp)
            dzdsc = dzdsc + 2.*real(q1*(0.,1.)*(j-1)*ctmp)
            dzdc2 = dzdc2 - 2.*real(q0*(j-1)**2*ctmp)
         endif
      enddo

C     compute jacobian and its first order derivatives
      jac    = (drds*dzdc - drdc*dzds)*r0
      djacds = drds2*dzdc*r0 + drds*dzdsc*r0 + dzdc*drds**2 
     &        -drdsc*dzds*r0 - drdc*dzds2*r0 - drdc*dzds*drds
      djacdc = drdsc*dzdc*r0 + drds*dzdc2*r0 + drds*dzdc*drdc
     &        -drdc2*dzds*r0 - drdc*dzdsc*r0 - dzds*drdc**2 

C     compute b10,b20,b30 and all first-order derivatices w.r.t.
C     (s,chi,phi) at (s0,chi0)
      if (s0.lt.1..and.i.gt.1) i=i-1 

C     h## for b10
      h1   = cs(i,n)    - s0
      h2   = cs(i+1,n)  - s0
      h3   = cs(i+2,n)  - s0
      d0   = (h1-h2)*(h2-h3)*(h3-h1)
      d1   = d0
      h01  = h3*h2*(h3-h2)/d0 
      h02  = h1*h3*(h1-h3)/d0
      h03  = h2*h1*(h2-h1)/d0
      h11  = (h2-h3)*(h2+h3)/d1
      h12  = (h3-h1)*(h3+h1)/d1
      h13  = (h1-h2)*(h1+h2)/d1

C     g## for b20 and b30
      h1   = csm(i,n)    - s0
      h2   = csm(i+1,n)  - s0
      h3   = csm(i+2,n)  - s0
      d0   = (h1-h2)*(h2-h3)*(h3-h1)
      d1   = d0
      g01  = h3*h2*(h3-h2)/d0 
      g02  = h1*h3*(h1-h3)/d0
      g03  = h2*h1*(h2-h1)/d0
      g11  = (h2-h3)*(h2+h3)/d1
      g12  = (h3-h1)*(h3+h1)/d1
      g13  = (h1-h2)*(h1+h2)/d1

      b10    = (0.,0.)
      b20    = (0.,0.)   
      b30    = (0.,0.)   
      db10ds = (0.,0.)
      db10dc = (0.,0.)
      db10dp = (0.,0.)
      db20ds = (0.,0.)
      db20dc = (0.,0.)
      db20dp = (0.,0.)
      db30ds = (0.,0.)
      db30dc = (0.,0.)
      db30dp = (0.,0.)

      do j=1,mmaxp(n)
         ctmp = exp((0.,1.)*(j-1+m1(n))*chi0)

         y1 = b1(i,j,n)
         y2 = b1(i+1,j,n)
         y3 = b1(i+2,j,n)
         y0 = h01*y1 + h02*y2 + h03*y3
         yp = h11*y1 + h12*y2 + h13*y3
         b10    = b10    + y0*ctmp
         db10ds = db10ds + yp*ctmp
         db10dc = db10dc + y0*ctmp*(0.,1.)*(j-1+m1(n))

         y1 = b2(i,j,n)
         y2 = b2(i+1,j,n)
         y3 = b2(i+2,j,n)
         y0 = g01*y1 + g02*y2 + g03*y3
         yp = g11*y1 + g12*y2 + g13*y3
         b20    = b20    + y0*ctmp
         db20ds = db20ds + yp*ctmp
         db20dc = db20dc + y0*ctmp*(0.,1.)*(j-1+m1(n))

         y1 = b3(i,j,n)
         y2 = b3(i+1,j,n)
         y3 = b3(i+2,j,n)
         y0 = g01*y1 + g02*y2 + g03*y3
         yp = g11*y1 + g12*y2 + g13*y3
         b30    = b30    + y0*ctmp
         db30ds = db30ds + yp*ctmp
         db30dc = db30dc + y0*ctmp*(0.,1.)*(j-1+m1(n))
      enddo
      db10dp =-(0.,1.)*nn(n)*b10
      db20dp =-(0.,1.)*nn(n)*b20
      db30dp =-(0.,1.)*nn(n)*b30
      
C     compute Br,phi,z-field at (r,phi,z)
C     and all the first-order derivatives w.r.t. (r,phi,z)
      br = (b10*drds+b20*drdc)/jac
      bp =-b30*r0/jac
      bz = (b10*dzds+b20*dzdc)/jac

      dbrds = db10ds*drds/jac + b10*drds2/jac - b10*drds*djacds/jac**2
     &       +db20ds*drdc/jac + b20*drdsc/jac - b20*drdc*djacds/jac**2 
      dbrdc = db10dc*drds/jac + b10*drdsc/jac - b10*drds*djacdc/jac**2
     &       +db20dc*drdc/jac + b20*drdc2/jac - b20*drdc*djacdc/jac**2 
      dbzds = db10ds*dzds/jac + b10*dzds2/jac - b10*dzds*djacds/jac**2
     &       +db20ds*dzdc/jac + b20*dzdsc/jac - b20*dzdc*djacds/jac**2 
      dbzdc = db10dc*dzds/jac + b10*dzdsc/jac - b10*dzds*djacdc/jac**2
     &       +db20dc*dzdc/jac + b20*dzdc2/jac - b20*dzdc*djacdc/jac**2 
      dbpds =-db30ds*r0/jac - b30*drds/jac + b30*r0*djacds/jac**2
      dbpdc =-db30dc*r0/jac - b30*drdc/jac + b30*r0*djacdc/jac**2

      dbrdr = (dbrds*dzdc-dbrdc*dzds)*r0/jac
      dbrdz = (dbrdc*drds-dbrds*drdc)*r0/jac
      dbrdp = (db10dp*drds+db20dp*drdc)/jac

      dbzdr = (dbzds*dzdc-dbzdc*dzds)*r0/jac
      dbzdz = (dbzdc*drds-dbzds*drdc)*r0/jac
      dbzdp = (db10dp*dzds+db20dp*dzdc)/jac

      dbpdr = (dbpds*dzdc-dbpdc*dzds)*r0/jac
      dbpdz = (dbpdc*drds-dbpds*drdc)*r0/jac
      dbpdp =-db30dp*r0/jac

C     compute brphz at real space point xrphz
      ctmp = 2.0*exp(-(0.,1.)*nn(n)*phi0)

      brphz(1) = brphz(1) + br*ctmp
      brphz(2) = brphz(2) + bp*ctmp
      brphz(3) = brphz(3) + bz*ctmp

      dbrphz(1,1) = dbrphz(1,1) + dbrdr*ctmp
      dbrphz(1,2) = dbrphz(1,2) + dbrdp*ctmp
      dbrphz(1,3) = dbrphz(1,3) + dbrdz*ctmp
      dbrphz(2,1) = dbrphz(2,1) + dbpdr*ctmp
      dbrphz(2,2) = dbrphz(2,2) + dbpdp*ctmp
      dbrphz(2,3) = dbrphz(2,3) + dbpdz*ctmp
      dbrphz(3,1) = dbrphz(3,1) + dbzdr*ctmp
      dbrphz(3,2) = dbrphz(3,2) + dbzdp*ctmp
      dbrphz(3,3) = dbrphz(3,3) + dbzdz*ctmp
      enddo
   
      return
      end

