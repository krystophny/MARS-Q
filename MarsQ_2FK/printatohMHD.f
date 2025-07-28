	SUBROUTINE printatoh(A,B,C,D,E,F,G,H,MX,MY,MR,NX,NY,NR9)
C
        USE DIMENSIM
	USE GLOBALM

	INTEGER MX,MY,MR,NX,NY,NR9
	COMPLEX*16 A(MX,MX,MR),B(MX,MX,MR),C(MX,MX,MR),D(MY,MY,MR),
     &             E(MX,MY,MR),F(MY,MX,MR),G(MY,MX,MR),H(MX,MY,MR)
        INTEGER IX,IY,IR
C
	LXCOL = (NX-1)*NXCOMP
	LYCOL = (NY-1)*NYCOMP
	
        if (1.eq.1) then
        IR=NR9
        write(*,*) 'A=['
        do IX=1,MX
          write(*,100) (A(IX,IY,IR),IY=1,MX)
        enddo
 100    format(800E16.8)

        write(*,*) ']; B=['
        do IX=1,MX
          write(*,100) (B(IX,IY,IR),IY=1,MX)
        enddo

        write(*,*) ']; C=['
        do IX=1,MX
          write(*,100) (C(IX,IY,IR),IY=1,MX)
        enddo

        write(*,*) ']; D=['
        do IX=1,MY
          write(*,100) (D(IX,IY,IR),IY=1,MY)
        enddo

        write(*,*) ']; E=['
        do IX=1,MX
          write(*,100) (E(IX,IY,IR),IY=1,MY)
        enddo

        write(*,*) ']; F=['
        do IX=1,MY
          write(*,100) (F(IX,IY,IR),IY=1,MX)
        enddo

        write(*,*) ']; G=['
        do IX=1,MY
          write(*,100) (G(IX,IY,IR),IY=1,MX)
        enddo

        write(*,*) ']; H=['
        do IX=1,MX
          write(*,100) (H(IX,IY,IR),IY=1,MY)
        enddo
        write(*,*) '];'

        else  

	WRITE(*,*) 'A...'
	DO IX=1,MX
	   WRITE(*,1) A(IX,KXB1+LXCOL,NR9),A(IX,KXV1+LXCOL,NR9),
     &                A(IX,KXJ2U+LXCOL,NR9),A(IX,KXJ3+LXCOL,NR9),
     &                A(IX,KXJ2L+LXCOL,NR9),A(IX,KXPD+LXCOL,NR9)
	ENDDO
 1	FORMAT(12E12.4)

	WRITE(*,*) 'B...'
	DO IX=1,MX
	   WRITE(*,1) B(IX,KXB1+LXCOL,NR9),B(IX,KXV1+LXCOL,NR9),
     &                B(IX,KXJ2U+LXCOL,NR9),B(IX,KXJ3+LXCOL,NR9),
     &                B(IX,KXJ2L+LXCOL,NR9),B(IX,KXPD+LXCOL,NR9)
	ENDDO

	WRITE(*,*) 'C...'
	DO IX=1,MX
	   WRITE(*,1) C(IX,KXB1+LXCOL,NR9),C(IX,KXV1+LXCOL,NR9),
     &                C(IX,KXJ2U+LXCOL,NR9),C(IX,KXJ3+LXCOL,NR9),
     &                C(IX,KXJ2L+LXCOL,NR9),C(IX,KXPD+LXCOL,NR9)
	ENDDO


	WRITE(*,*) 'D...'
	DO IY=1,MY
	   WRITE(*,2) D(IY,KYV2+LYCOL,NR9),D(IY,KYV3+LYCOL,NR9),
     &                D(IY,KYJ1+LYCOL,NR9),D(IY,KYB3+LYCOL,NR9),
     &                D(IY,KYB2+LYCOL,NR9),D(IY,KYPR+LYCOL,NR9)
	ENDDO
 2	FORMAT(12E12.4)


	WRITE(*,*) 'E...'
	DO IX=1,MX 
	   WRITE(*,2) E(IX,KYV2+LYCOL,NR9),E(IX,KYV3+LYCOL,NR9),
     &                E(IX,KYJ1+LYCOL,NR9),E(IX,KYB3+LYCOL,NR9),
     &                E(IX,KYB2+LYCOL,NR9),E(IX,KYPR+LYCOL,NR9)
	ENDDO

	WRITE(*,*) 'F...'
	DO IY=1,MY
	   WRITE(*,1) F(IY,KXB1+LXCOL,NR9),F(IY,KXV1+LXCOL,NR9),
     &                F(IY,KXJ2U+LXCOL,NR9),F(IY,KXJ3+LXCOL,NR9),
     &                F(IY,KXJ2L+LXCOL,NR9),F(IY,KXPD+LXCOL,NR9)
	ENDDO

	WRITE(*,*) 'G...'
	DO IY=1,MY
	   WRITE(*,1) G(IY,KXB1+LXCOL,NR9),G(IY,KXV1+LXCOL,NR9),
     &                G(IY,KXJ2U+LXCOL,NR9),G(IY,KXJ3+LXCOL,NR9),
     &                G(IY,KXJ2L+LXCOL,NR9),G(IY,KXPD+LXCOL,NR9)
	ENDDO

	WRITE(*,*) 'H...'
	DO IX=1,MX 
	   WRITE(*,2) H(IX,KYV2+LYCOL,NR9),H(IX,KYV3+LYCOL,NR9),
     &                H(IX,KYJ1+LYCOL,NR9),H(IX,KYB3+LYCOL,NR9),
     &                H(IX,KYB2+LYCOL,NR9),H(IX,KYPR+LYCOL,NR9)
	ENDDO

C	WRITE(*,'("NX,NY,NR=",3I5)')NX,NY,NR
C       WRITE(*,'("A,B     =",1P4E12.5)')A(NX,NY,NR),B(NX,NY,NR)
C       WRITE(*,'("C,D     =",1P4E12.5)')C(NX,NY,NR),D(NX,NY,NR)
C       WRITE(*,'("E,F     =",1P4E12.5)')E(NX,NY,NR),F(NX,NY,NR)
C       WRITE(*,'("G,H     =",1P4E12.5)')G(NX,NY,NR),H(NX,NY,NR)
C
        endif

	RETURN
        END
