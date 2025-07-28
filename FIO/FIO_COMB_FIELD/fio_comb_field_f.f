C     ========================================================
C     read field data from different rows of coils and combine 
C     them into a single file, using given coil phasing, 
C     YQ Liu, 11/06/2017
C     =======================================================
      program fio_comb_field_f

      implicit none
      integer::i,j,k,n,ns,m1,m2
      integer::nsp,nsv,mmaxe,mmaxp
      real*8::pi_value,pu,pm,pl,fu,fm,fl
      complex*16::epu,epm,epl,bc
      real*8,dimension(8)::rtmp
      real*8,dimension(6)::ru,rm,rl
      character(len=128) dn,fn1,fn2,fn3,fn4
      
      dn  = '/home/liuy/Work/ITER/Data/'
C     fn1 = 'BPLASMA_MARSF_n6_cU_Pr03_tfte065.IN'
C     fn2 = 'BPLASMA_MARSF_n6_cM_Pr03_tfte065.IN'
C     fn3 = 'BPLASMA_MARSF_n6_cL_Pr03_tfte065.IN'
C     fn4 = 'BPLASMA_MARSF_n6_cC_Pr03_tfte065.IN'
      fn1 = 'BPLASMA_MARSF_n3_cU_VAC.IN'
      fn2 = 'BPLASMA_MARSF_n3_cM_VAC.IN'
      fn3 = 'BPLASMA_MARSF_n3_cL_VAC.IN'
      fn4 = 'BPLASMA_MARSF_n3_cC_VAC.IN'

      fn1 = trim(dn)//trim(fn1)
      fn2 = trim(dn)//trim(fn2)
      fn3 = trim(dn)//trim(fn3)
      fn4 = trim(dn)//trim(fn4)

      pi_value = acos(-1.)

      pu  = pi_value*111.9/180.0
      pm  = pi_value*0.0/180.0
      pl  = pi_value*267.9/180.0

      epu = exp((0.,1.)*pu) 
      epm = exp((0.,1.)*pm) 
      epl = exp((0.,1.)*pl) 

C     open files to read/write field data
      open(11,file=trim(fn1),status='old',form='formatted')
      open(12,file=trim(fn2),status='old',form='formatted')
      open(13,file=trim(fn3),status='old',form='formatted')
      open(14,file=trim(fn4),form='formatted')

C     read/write dimensions       
      read(11,*)    n,mmaxe,m1,m2,nsp,nsv,i,rtmp(1),rtmp(2)
      read(12,*)    n,mmaxe,m1,m2,nsp,nsv,i,rtmp(1),rtmp(2)
      read(13,*)    n,mmaxe,m1,m2,nsp,nsv,i,rtmp(1),rtmp(2)
      write(14,118) n,mmaxe,m1,m2,nsp,nsv,i,rtmp(1),rtmp(2)
 118  format(7(I4,1X),2(E17.10,1X))
 119  format(8(E17.10,1X))

      if (abs(n).eq.3) then
         fu = 0.6645
         fm = 0.4968
         fl = 0.6840
      elseif (abs(n).eq.4) then
         fu = 0.6126
         fm = 0.4774
         fl = 0.6264
      elseif (abs(n).eq.5) then
         fu = 0.5494
         fm = 0.4530
         fl = 0.5565
      elseif (abs(n).eq.6) then
         fu = 0.4772
         fm = 0.4243
         fl = 0.4773
      endif

      mmaxp = m2-m1 + 1
      ns    = nsp + nsv

C     read/write radial mesh and q-profile data 
      do i=1,ns
         read(11,*)    rtmp(1),rtmp(2),rtmp(3)
         read(12,*)    rtmp(1),rtmp(2),rtmp(3)
         read(13,*)    rtmp(1),rtmp(2),rtmp(3)
         write(14,119) rtmp(1),rtmp(2),rtmp(3)
      enddo

C     read/write coordinates mapping data
      do j=1,mmaxe
      do i=1,ns
         read(11,*)    (rtmp(k),k=1,8)
         read(12,*)    (rtmp(k),k=1,8)
         read(13,*)    (rtmp(k),k=1,8)
         write(14,119) (rtmp(k),k=1,8)
      enddo
      enddo

C     read/write b-field data
      do j=1,mmaxp
      do i=1,ns
         read(11,*) (ru(k),k=1,6)
         read(12,*) (rm(k),k=1,6)
         read(13,*) (rl(k),k=1,6)
         bc = cmplx(ru(1),ru(2))*epu*fu + cmplx(rm(1),rm(2))*epm*fm +
     &        cmplx(rl(1),rl(2))*epl*fl
         rtmp(1) = real(bc)
         rtmp(2) = imag(bc)
         bc = cmplx(ru(3),ru(4))*epu*fu + cmplx(rm(3),rm(4))*epm*fm +
     &        cmplx(rl(3),rl(4))*epl*fl
         rtmp(3) = real(bc)
         rtmp(4) = imag(bc)
         bc = cmplx(ru(5),ru(6))*epu*fu + cmplx(rm(5),rm(6))*epm*fm +
     &        cmplx(rl(5),rl(6))*epl*fl
         rtmp(5) = real(bc)
         rtmp(6) = imag(bc)
         write(14,119) (rtmp(k),k=1,6)
      enddo
      enddo

C     read/write x-field data
      do j=1,mmaxp
      do i=1,nsp+1
         read(11,*) (ru(k),k=1,2)
         read(12,*) (rm(k),k=1,2)
         read(13,*) (rl(k),k=1,2)
         bc = cmplx(ru(1),ru(2))*epu*fu + cmplx(rm(1),rm(2))*fm +
     &        cmplx(rl(1),rl(2))*epl*fl
         write(14,119) real(bc),imag(bc)
      enddo
      enddo

      close(11)
      close(12)
      close(13)
      close(14)

      return 
      end

