      module fio_eval_field_fm

      integer::nr,nz,ns,mmaxe,nmax
      integer,dimension(:),allocatable::nn,mmaxp,m1,m2
      real*8 ::raxis,zaxis,pi_value
      real*8,dimension(:),  allocatable::cs,csm
      real*8,dimension(:,:),allocatable::r,z,s,chi
      complex*16,dimension(:,:),allocatable::rmi,zmi,rmm,zmm
      complex*16,dimension(:,:,:),allocatable::b1,b2,b3,x1

      end module fio_eval_field_fm
