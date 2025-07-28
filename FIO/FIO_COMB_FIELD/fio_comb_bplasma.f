C     ========================================================
C     read field data from BPLASMA_MARSF_* file and 
C     output into BPLASMA.OUT format
C     YQ Liu, 11/13/2017
C     =======================================================
      program fio_comb_bplasma

      implicit none
      integer::i,j,k,n,ns,m1,m2
      integer::nsp,nsv,mmaxp,mmaxe
      real*8,dimension(8)::rtmp
      character(len=128) dn,fn1,fn2
      
      dn  = '/home/ITER/liuy2/ITERout/case5/'
C     fn1 = 'BPLASMA_MARSF_n3_cC_Pr03_tfte2.IN'
C     fn2 = 'BPLASMA_n3_cC_Pr03_tfte2.OUT'
      fn1 = 'BPLASMA_MARSF_n3_cC_vac.IN'
      fn2 = 'BPLASMA_n3_cC_vac.OUT'

      fn1 = trim(dn)//trim(fn1)
      fn2 = trim(dn)//trim(fn2)

C     open files to read/write field data
      open(11,file=trim(fn1),status='old',form='formatted')
      open(12,file=trim(fn2),form='formatted')

C     read/write dimensions       
      read(11,*)    n,mmaxe,m1,m2,nsp,nsv,i,rtmp(1),rtmp(2)

      mmaxp = m2-m1 + 1
      ns    = nsp + nsv

      write(12,1171) mmaxp,ns,dfloat(n),0,0,0
      do j=1,mmaxp
         write(12,1172) dfloat(j-1+m1),dfloat(j-1+m1),dfloat(j-1+m1),
     &                  dfloat(j-1+m1),dfloat(j-1+m1),dfloat(j-1+m1)
      enddo   

 1171 FORMAT(I5,1X,I5,1X,E9.2,3(1X,I2))
 1172 FORMAT(E16.8E3,5(1X,E16.8E3))

C     read radial mesh and q-profile data 
      do i=1,ns
         read(11,*)    rtmp(1),rtmp(2),rtmp(3)
      enddo

C     read coordinates mapping data
      do j=1,mmaxe
      do i=1,ns
         read(11,*)    (rtmp(k),k=1,8)
      enddo
      enddo

C     read/write b-field data
      do j=1,mmaxp
      do i=1,ns
         read(11,*)     (rtmp(k),k=1,6)
         write(12,1172) (rtmp(k),k=1,6)
      enddo
      enddo

      close(11)
      close(12)

      return 
      end

