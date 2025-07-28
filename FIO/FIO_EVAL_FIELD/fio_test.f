      program field_test

      implicit none
      integer::B_FIELD,A_FIELD,X_FIELD,ierr
      real*8,dimension(3)::x,y

C     define field type
      B_FIELD = 1
      A_FIELD = 2
      X_FIELD = 3

C     read in MARS-F field data on a grid
      call fio_init_field_f(ierr)

C     specify point x in (r,phi,z) space
      x(1) = 0.69656250E+00   !r [m]
      x(2) = 0.    !phi [rad]
      x(3) = 0.17187500E-01   !z [m]

C     compute y=B-field at point x
      call fio_eval_field_f(B_FIELD,x,y,ierr)

      write(*,'("xrphz = ",3(e15.8,1x))') x(1),x(2),x(3)
      write(*,'("brphz = ",3(e15.8,1x))') y(1),y(2),y(3)

C     free allocated memory
      call fio_finish_field_f

      return 
      end

