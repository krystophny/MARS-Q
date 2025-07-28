      module fio_eval_field_fm

      integer::nr,nz,ns,n,mmaxe,mmaxp,m1,m2
      real*8 ::raxis,zaxis
      real*8,dimension(:),  allocatable::cs,csm
      real*8,dimension(:,:),allocatable::r,z,s,chi
      complex*16,dimension(:,:),allocatable::rmi,zmi,rmm,zmm,b1,b2

      end module fio_eval_field_fm
