      program field_test

      implicit none
      integer::B_FIELD,A_FIELD,X_FIELD,ierr
      real*8,dimension(3)::x
      complex*16,dimension(3)::y

C     define field type
      B_FIELD = 1
      A_FIELD = 2
      X_FIELD = 3

C     read in MARS-F field data on a grid
      call fio_init_field_f(ierr)

C     specify point x in (r,phi,z) space
      x(1) = 6.2273   !r [m]
      x(2) = 0.    !phi [rad]
      x(3) = 0.0571   !z [m]

C     compute y=B-field at point x
      call fio_eval_field_f(B_FIELD,x,y,ierr)

      write(*,'("xrphz = ",3(e15.8,1x))') x(1),x(2),x(3)
      write(*,'("brphz = ",6(e15.8,1x))') real(y(1)),imag(y(1)),
     &                                    real(y(2)),imag(y(2)),
     &                                    real(y(3)),imag(y(3))

C     free allocated memory
      call fio_finish_field_f

      return 
      end

