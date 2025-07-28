      program field_test

      implicit none
      integer::B_FIELD,A_FIELD,ierr,k,i,j
      real*8 ::dr,dphi,dz,divb,modb
      real*8,dimension(3,7)::x,y
      real*8::rmin,rmax,zmin,zmax,nrr,nzz,hhr,hhz

C     define field type
      B_FIELD = 1
      A_FIELD = 2

C     specify mesh size for finite differencing
      dr   = 1.0e-4;
      dphi = 1.0e-4;   ! in [rad]
      dz   = 1.0e-4;
 
      call fio_init_field_f(ierr)

      
      rmin = 0.85
      rmax = 2.45
      zmin =-1.6
      zmax = 1.6
      nrr  = 110
      nzz  = 110

      hhr  = (rmax-rmin)/(nrr-1)
      hhz  = (zmax-zmin)/(nzz-1)

      open(99,file='BPLASMA_RECTRZ_FIO.OUT')
      rewind(99)

      do i=3,nrr-2
      do j=3,nzz-2

C      do i=1,1
C      do j=1,1

C     specify (R,PHI,Z) coordinates for 6 points (2-7) near the central point (1)
      x(1,1) = rmin + (i-1)*hhr   ! R-value of the central point
      x(2,1) = 1.0                ! PHI-value of the central point
      x(3,1) = zmin + (j-1)*hhz   ! Z-value of the central point

C      x(1,1) = 1.7530612
C      x(3,1) =-0.13877551

C      x(1,1) = 1.7530612
C      x(3,1) =-0.20408163

      do k=2,7
         x(:,k) = x(:,1);
      enddo

      x(1,2) = x(1,1) - dr
      x(1,3) = x(1,1) + dr

      x(2,6) = x(2,1) - dphi
      x(2,7) = x(2,1) + dphi

      x(3,4) = x(3,1) - dz
      x(3,5) = x(3,1) + dz

C     compute B-field values y at points x
      
      do k=1,7
         call fio_eval_field_f(B_FIELD,x(:,k),y(:,k),ierr)
      enddo


C     compute divB at the central point x(:,1)
C     and normalize by |B|=modb
      modb = sqrt(y(1,1)**2+y(2,1)**2+y(3,1)**2)
      divb = (x(1,3)*y(1,3)-x(1,2)*y(1,2))/2./dr/x(1,1) +
     &       (y(2,7)-y(2,6))/2./dphi/x(1,1) +
     &       (y(3,5)-y(3,4))/2./dz

C     write(*,*) 'test divB: rel.error=',divb/modb
      write(99,100) x(1,1),x(3,1),y(1,1),y(2,1),y(3,1),divb/modb
 100  format(6(e15.8,1x))

      enddo
      enddo

      close(99)

      call fio_finish_field_f

      return 
      end

