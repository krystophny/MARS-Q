      program fio_main

      implicit none
      integer::B_FIELD,k,ierr
      real*8,dimension(3)::x
      complex*16,dimension(3)::y
      complex*16,dimension(3,3)::yp

C     define field type
      B_FIELD = 1

C     read in MARS-F field data on a grid
      call fio_init_field_f(ierr)

C     specify point x in (r,phi,z) space
      x(1) = 3.00   !r [m]
      x(2) = 0.    !phi [rad]
      x(3) = 0.99  !z [m]

C     compute y=B-field at point x
      call fio_eval_field_f(B_FIELD,x,y,yp,ierr)

      write(*,'("xrphz = ",3(e15.8,1x))') x(1),x(2),x(3)
      write(*,'("brphz = ",3(e15.8,1x))') y(1),y(2),y(3)
      do k=1,3
         write(*,'("dbrphz = ",3(e15.8,1x))') yp(k,1),yp(k,2),yp(k,3)
      enddo

C     free allocated memory
      call fio_finish_field_f

      return 
      end

