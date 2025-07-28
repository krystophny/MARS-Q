      program field_test
      use fio_eval_field_fm

      implicit none
      integer::B_FIELD,A_FIELD,X_FIELD,ierr,i,j,k,NP
      real*8,dimension(3)::x,y
      real*8::rmin,rmax,zmin,zmax,nrr,nzz,hhr,hhz
      real*8::s0,chi0,phi0
      complex*16::ctmp

C     define field type
      B_FIELD = 1
      A_FIELD = 2
      X_FIELD = 3

C     read in MARS-F field data on a grid
      call fio_init_field_f(ierr)

C     define rectangular mesh in (R,Z)      
      rmin = 3.52
      rmax = 8.88
      zmin =-5.28
      zmax = 5.28
      nrr  = 101
      nzz  = 201

      hhr  = (rmax-rmin)/(nrr-1)
      hhz  = (zmax-zmin)/(nzz-1)

C     compute normal displacement at rectangular mesh
      open(99,file='XPLASMA_RECTRZ_FIO.OUT')
      rewind(99)

      do i=1,nrr
      do j=1,nzz
         x(1) = rmin + (i-1)*hhr   
         x(2) = 0.0                
         x(3) = zmin + (j-1)*hhz

         call fio_eval_field_f(X_FIELD,x,y,ierr)
         write(99,100) x(1),x(3),y(1)
100      format(3(e15.8,1x))
      enddo   
      enddo   

      close(99)

C     compute normal displacement along surface s=s0
      open(99,file='XPLASMA_SURFACE_FIO.OUT')
      rewind(99)

      s0   = 0.995
      phi0 = 0.
      NP   = 201

C     search for surface number k: cs(k)=s0
      k  = 0
      do i=1,ns-1
         if ( (cs(i)-s0)*(cs(i+1)-s0).le.0. ) k=i
      enddo
      if (k.eq.0) then
         ierr = 3
         write(*,*) 'fio_test3: ierr=', ierr
         write(*,*) 's0=',s0
         stop
      endif
      i = k

      do k=1,NP-1
      chi0 = -pi_value + (k-1)*2*pi_value/float(NP-1)
      x(1) = 0.
      x(2) = phi0
      x(3) = 0.
      do j=1,mmaxe
         if (j.eq.1) then
            x(1) = x(1) + real(rmm(i,j))
            x(3) = x(3) + real(zmm(i,j))
         else
            ctmp = exp((0.,1.)*(j-1)*chi0)
            x(1) = x(1) + 2.*real(rmm(i,j)*ctmp)
            x(3) = x(3) + 2.*real(zmm(i,j)*ctmp)
         endif
      enddo   

      call fio_eval_field_f(X_FIELD,x,y,ierr)
      write(99,200) chi0,x(1),x(2),x(3),y(1)
200   format(5(e15.8,1x))
      enddo   

      close(99)

C     free allocated memory
      call fio_finish_field_f

      return 
      end

