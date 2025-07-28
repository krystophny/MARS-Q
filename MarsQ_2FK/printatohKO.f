	SUBROUTINE printatoh(A,B,C,D,E,F,G,H,MX,MY,MR,NX,NY,NR9)
C
        USE DIMENSIM
	USE GLOBALM

	INTEGER MX,MY,MR,NX,NY,NR9
	COMPLEX*16 A(MX,MX,MR),B(MX,MX,MR),C(MX,MX,MR),D(MY,MY,MR),
     &             E(MX,MY,MR),F(MY,MX,MR),G(MY,MX,MR),H(MX,MY,MR)
        INTEGER IX,IY,IR,IZ
C
	LYCOL = (NY-1)*NYCOMP
	
	WRITE(*,*) 'D...'
	DO IZ=1,MSMAX
           IY = (IZ-1)*NYCOMP
	   WRITE(*,2) D(KYPPARA+LYCOL,KYX2+IY,NR9),
     &                D(KYPPARA+LYCOL,KYB2+IY,NR9),
     &                D(KYPPARA+LYCOL,KYB3+IY,NR9),
     &                D(KYPPERP+LYCOL,KYX2+IY,NR9),
     &                D(KYPPERP+LYCOL,KYB2+IY,NR9),
     &                D(KYPPERP+LYCOL,KYB3+IY,NR9)
	ENDDO
 2	FORMAT(12E12.4)


	WRITE(*,*) 'F...'
	DO IZ=1,MSMAX
           IX = (IZ-1)*NXCOMP
	   WRITE(*,1) F(KYPPARA+LYCOL,KXX1+IX,NR9),
     &                F(KYPPARA+LYCOL,KXB1+IX,NR9),
     &                F(KYPPERP+LYCOL,KXX1+IX,NR9),
     &                F(KYPPERP+LYCOL,KXB1+IX,NR9)
	ENDDO
 1	FORMAT(8E12.4)

	WRITE(*,*) 'G...'
	DO IZ=1,MSMAX
           IX = (IZ-1)*NXCOMP
	   WRITE(*,1) G(KYPPARA+LYCOL,KXX1+IX,NR9),
     &                G(KYPPARA+LYCOL,KXB1+IX,NR9),
     &                G(KYPPERP+LYCOL,KXX1+IX,NR9),
     &                G(KYPPERP+LYCOL,KXB1+IX,NR9)
	ENDDO

	RETURN
        END
