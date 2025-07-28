      program test_NTV 
C**********************************************
C*     test the NTV torque calculation
C*
C*      by Y. Sun
C**********************************************
!     use NTV
!     use NTV_lib

!     call test1
!     call test2
      call test3
      end

      subroutine test3
C***************************************************
C    test the NTV modeling
C    use the data from MARS-F output
C***************************************************
      use NTV           ! for NTV torque calculation
      use common_ntv    ! for common variants needed for NTV calculation
      integer i

      call basic_input    ! set dimensions of the arrays
      call ntv_predata    ! load data from MARS-F output
      call gntv1d         ! calculate NTV torque density
! ------ save results
        open(21,file='temp/spec.dat') 
        open(22,file='temp/prof.dat') 
 
         do i=1,Nr
          write(22,99)nu_prof(i,:)
         end do
         close(22)

        write(21,99)rho
         do i=1,(2*m_max+1)
          write(21,99)sqrt(bc(i,1,:)**2+bs(i,1,:)**2)
         end do
         close(21)

	open(31,file='temp/test3.dat') 
         do i=1,nr
          write(31,99)rho(i),T_ntv(i)
         end do
         close(31)

99     format(1x,500E15.7)

      end  subroutine test3

      subroutine test2
C***************************************************
C    test the pitch angle integral
C      energy(x) dependence of I_kn.
C     for 1 collsionality case
C***************************************************
      use NTV
	use NTV_lib
	implicit none
        integer,parameter::nk=101,nt=201,nx=81,ic=1,n_m=3
        integer mm(n_m),i,n
	real bc(n_m),bs(n_m)
        real q,qwe,wB0,qws,qwsT,wti,rho_i,R0,eps,norm_Z
        real nud0,x(nx),nud(nx),I_kn(nx),lambda(3)
        real T_ntv
        integer,parameter:: n_max=500
        integer n_k02
        real x_k02(n_max),k02(n_max)
     
        open(21,file='temp/test2.dat') 
        call test_spec(n,n_m,mm,bc,bs)
        print*,'n=',n
        print*,'mm=',mm
	print*,'bc=',bc
	print*,'bs=',bs
        call test_ex2(n,nx,x,nud,nud0,
     .    q,qwe,wB0,qws,qwsT,wti,rho_i,R0,eps,norm_Z)
        print*,'nu*d=',nud0
        print*,'x=',x
        print*,'nud=',nud
	print*,'q=',q,'wB0=',wB0,'qwe=',qwe

!------- load data for the calculation of resonant pitch
        call cal_resk(n_k02,x_k02,k02)
          if (n_max.LT.n_k02)then
             stop 'dimensional error, (n_max<n_k02) '
          end if

      call gntv(n_m,mm,bc,bs,q,n,nk,nx,n_k02,x_k02(1:n_k02)
     .  ,k02(1:n_k02),x,nud,qwe,wB0,qws,qwsT,wti,rho_i,R0,eps,norm_Z
     .                   ,ic,I_kn,lambda,T_ntv) 
       print*,'      nud0','      I_kn'
             do i=1,nx
                print*,i,x(i),I_kn(i)
               write(21,99)x(i),I_kn(i)
              end do
              print*,'lambda',lambda

	close(21)

99     format(1x,5E15.7)
      end subroutine test2

      subroutine test1
C***************************************************
C    test the pitch angle integral
C        collisionality dependence of I_kn
C        for the same energy of the particles
C***************************************************
       use NTV
	use NTV_lib
	implicit none
        integer,parameter::nk=101,nt=201,nx=101,ic=1,n_m=3
        integer mm(n_m),i,n
	real bc(n_m),bs(n_m)
        real q,qwe,wB0,qws,qwsT,wti,rho_i,R0,eps,norm_Z
        real nud0(nx),x(nx),nud(nx),I_kn(nx),lambda(3)
        real T_ntv
        integer,parameter:: n_max=500
        integer n_k02
        real x_k02(n_max),k02(n_max)
     
        open(21,file='temp/test1.dat') 
        call test_spec(n,n_m,mm,bc,bs)

        print*,'n=',n
        print*,'mm=',mm
	print*,'bc=',bc
	print*,'bs=',bs

        call test_ex1(n,nx,x,nud,nud0,
     .    q,qwe,wB0,qws,qwsT,wti,rho_i,R0,eps,norm_Z)
        print*,'nu*d=',nud0
        print*,'x=',x
        print*,'nud=',nud
	print*,'q=',q,'wB0=',wB0,'qwe=',qwe

!------- load data for the calculation of resonant pitch
        call cal_resk(n_k02,x_k02,k02)
          if (n_max.LT.n_k02)then
             stop 'dimensional error, (n_max<n_k02) '
          end if

      call gntv(n_m,mm,bc,bs,q,n,nk,nx,n_k02,x_k02(1:n_k02)
     .  ,k02(1:n_k02),x,nud,qwe,wB0,qws,qwsT,wti,rho_i,R0,eps,norm_Z
     .                   ,ic,I_kn,lambda,T_ntv) 
       print*,'      nud0','      I_kn'
             do i=1,nx
                print*,nud0(i),I_kn(i)
               write(21,99)nud0(i),I_kn(i)
              end do
              print*,'lambda',lambda

	close(21)

99     format(1x,5E15.7)
      end subroutine test1

      subroutine cal_resk(n_k02,x_k02,k02)
C***************************************************
C  load data for the calculation of resonant pitch
C***************************************************
      implicit none
      real x_k02(*),k02(*)
      integer n_k02,i
      real x,k
!     real,allocatable::t_x(:),t_k(:)

      open(11,file='temp/kappa0.dat',status='old') 
      read(11,*)n_k02
      do i=1,n_k02
         read(11,*)x_k02(i),k02(i)
      end do
      close(11)
      end subroutine cal_resk

      subroutine test_spec(n,n_m,mm,bc,bs)
C***************************************************
C    test spec profiles
C***************************************************
      integer n,n_m,mm(n_m)
      real bc(n_m),bs(n_m)
        n=1
	mm(1)=1
	mm(2)=2
        mm(3)=3
	bc=0.0
	bs=0.0
	bc(1)=1.0e-3
	bs(2)=0.0e-3
      end subroutine test_spec

      subroutine test_ex1(n,nx,x,nud,nud0,
     .    q,qwe,wB0,qws,qwsT,wti,rho_i,R0,eps,norm_Z)
C***************************************************
C    test plasma profiles
C***************************************************
       use NTV
	implicit none
       integer nx,Z,n
       real q,qwe,wB0,qws,qwsT,wti,rho_i,R0,eps,norm_Z
       real nud0(nx),x(nx),k_ind(nx),c0,f_nu(nx),nud(nx)
       real nui(nx)
       Z=1
        q=1.1
        qwe=2.0e+3
        wB0=sign(0.5e-0*abs(qwe),1.0*Z)

        qws=0.0
        qwsT=0.0
        wti=0.0
        rho_i=1.0
        R0=1.0
        eps=0.1
        norm_Z=1.0
        x=1.0
	call xgrid_f(nx,-3.0,2.0,k_ind)
        nud0=10.0**k_ind

        if (Z.GT.0)then
          c0=1.0
        else
          c0=sqrt(1836.0)
        end if
        nui=2.0*nud0*abs(n*qwe)*eps*c0
        call f_nud(x,nx,Z,f_nu)
        nud=0.5*nui/eps*f_nu 

       end subroutine test_ex1

      subroutine test_ex2(n,nx,x,nud,nud0,
     .    q,qwe,wB0,qws,qwsT,wti,rho_i,R0,eps,norm_Z)
C***************************************************
C    test plasma profiles
C***************************************************
       use NTV
	implicit none
       integer nx,Z,n,i
       real q,qwe,wB0,qws,qwsT,wti,rho_i,R0,eps,norm_Z
       real nud0,x(nx),k_ind(nx),c0,f_nu(nx),nud(nx)
       real nui
        Z=1
        q=1.1
        qwe=2.0e+3
        wB0=sign(0.5e-0*abs(qwe),1.0*Z)

        qws=0.0
        qwsT=0.0
        wti=0.0
        rho_i=1.0
        R0=1.0
        eps=0.1
        norm_Z=1.0
        call xgrid_f(nx,1.0e-2,2.0e1,x)

        nud0=1.0e-2

        if (Z.GT.0)then
          c0=1.0
        else
          c0=sqrt(1836.0)
        end if
        nui=2.0*nud0*abs(n*qwe)*eps*c0
        call f_nud(x,nx,Z,f_nu)
        nud=0.5*nui/eps*f_nu 

       end subroutine test_ex2
