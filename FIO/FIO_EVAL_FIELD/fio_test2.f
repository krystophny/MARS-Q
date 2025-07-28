      program field_test

      implicit none
      integer::B_FIELD,A_FIELD,X_FIELD,ierr
      real*8,dimension(3)::x,y
      integer::NR,NZ,j,k

C     define field type
      B_FIELD = 1
      A_FIELD = 2
      X_FIELD = 3

C     read in MARS-F field data on a grid
      call fio_init_field_f(ierr)

C     specify point x in (r,phi,z) space
      x(1) = 1.8   !r [m]
      x(2) = 0.    !phi [rad]
      x(3) = 0.    !z [m]

      NR = 101
      NZ = 101

      open(10,file='BFIELD.OUT',status='new')
      rewind(10)
 
      do j=1,NR
      do k=1,NZ
         x(1) = 1.0 + 1.4*(j-1)/(NR-1)
         x(2) = 0.
         x(3) = -1.2 + 1.0*(k-1)/(NZ-1)

C     compute y=B-field at point x
      call fio_eval_field_f(B_FIELD,x,y,ierr)

C     write(*,'("xrphz = ",3(e15.8,1x))') x(1),x(2),x(3)
C     write(*,'("brphz = ",3(e15.8,1x))') y(1),y(2),y(3)

         write(10,119) x(1),x(3),y(1),y(2),y(3)
      enddo
      enddo
 119  FORMAT(5(E12.5,1X))
      
      close(10)

C     free allocated memory
      call fio_finish_field_f

      return 
      end

