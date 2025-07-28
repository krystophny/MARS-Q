c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c--- OUTPUT EQUILIBRIUM QUANTITIES FOR T7 --- Y.Q. Liu 13.09.2012 ------
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

      Subroutine OutputT7           
      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      IMPLICIT NONE
      INCLUDE 'comioc.inc'
      INCLUDE 'comfft.inc'

      INTEGER I,J,M,IQ0
      INTEGER NPSTRT,MOUTP1
C
C     DECLARE QUANTITIES TO STORE OUTPUT FOR T7
      REAL*8 T7TMPS,T7EPS0,T7D2QDR20,T7D2PDR20,CS00,
     &       DQDU0,D2QDU20,DPDU0,D2PDU20,
     &       T7SA,T7SB,T7SC,T7SD,T7SE
      REAL*8,DIMENSION(:),ALLOCATABLE::T7R,T7DUDR,T7DQDR,T7DPDR,T7F,
     &       T7DGDR,T7DFRDR, T7D2GDR2,T7D2GC0DR2,T7D2PD0DR2,T7DC0DR,
     &       T7TMPA,T7TMPB
      REAL*8,DIMENSION(:,:),ALLOCATABLE::WORKS,RW1,RW2,RW3,
     &       T7AK,T7BK,T7CK,T7DK,T7EK,T7FK,T7GK,T7DADRK
      COMPLEX*16,DIMENSION(:,:),ALLOCATABLE::CW1

      ALLOCATE(
     &   T7R(NRP1),    T7TMPA(NRP1),    T7DUDR(NRP1),  T7DQDR(NRP1),
     &   T7TMPB(NRP1), T7DPDR(NRP1),    T7F(NRP1),     T7DGDR(NRP1),
     &   T7DFRDR(NRP1),T7D2GDR2(NRP1),  T7D2GC0DR2(NRP1), 
     &   T7D2PD0DR2(NRP1),    T7DC0DR(NRP1),
     &   RW1(NRP1,NCHI),      RW2(NRP1,NCHI),      RW3(NRP1,NCHI),
     &   WORKS(NRP1,2*NCHI),  CW1(NRP1,MOUTPUT+1),
     &   T7AK(NRP1,0:MOUTPUT),T7BK(NRP1,0:MOUTPUT),T7CK(NRP1,0:MOUTPUT),
     &   T7DK(NRP1,0:MOUTPUT),T7EK(NRP1,0:MOUTPUT),T7FK(NRP1,0:MOUTPUT),
     &   T7GK(NRP1,0:MOUTPUT),T7DADRK(NRP1,0:MOUTPUT)
     &        )

      INCLUDE 'setfft.inc'
      SUBNAM    = 'OutputT7'
C
      NPSTRT    =  1
      MOUTP1    = MOUTPUT+1
C
      CS00   = CS(1)
      CS(1)  = 0.
      CSH(1) = CS(2) - CS(1)

C     COMPUTE MINOR RADIUS T7R
      T7R(1) = 0.
      DO I=2,NRP1
         T7R(I) = T7R(I-1)+2.*CSM(I-1)*QM(I-1)/TM(I-1)*CSH(I-1)
      ENDDO
      T7TMPS = T7R(NRP1)
      DO I=2,NRP1
         T7R(I) = T7R(I)/T7TMPS
      ENDDO

C     COMPUTE DS/DR AT RATIONAL SURFACES, WHERE R=SQRT(T7R)
      WRITE(*,*) 'OutputT7: DSDR AT RATIONAL SURFACES:'
      DO I=1,NRATSURF
         J=IRATSURF(I)
         WRITE(*,110) J,CS(J),T7TMPS*SQRT(T7R(J))*T(J)/CS(J)/Q(J)
      ENDDO
 110  FORMAT(I4,2(1X,E13.5))

C     COMPUTE T7A=a/R
      T7EPS0 = SQRT(DPSIDS(NRP1)*T7TMPS)

C     COMPUTE Q(1),DQDU0,D2QDU20, WHERE U=S^2
C     USING THREE POINTS Q(2),Q(3) AND Q(4)
      IQ0  = 6
      T7SA = CS(IQ0+2)**2 - CS(IQ0+1)**2
      T7SB = .5*(CS(IQ0+2)**4-CS(IQ0+1)**4)
      T7SC = CS(IQ0+3)**2 - CS(IQ0+1)**2
      T7SD = .5*(CS(IQ0+3)**4-CS(IQ0+1)**4)
      T7SE = T7SA*T7SD - T7SC*T7SB
      DQDU0   = ((Q(IQ0+2)-Q(IQ0+1))*T7SD-(Q(IQ0+3)-Q(IQ0+1))*T7SB)/T7SE
      D2QDU20 = ((Q(IQ0+3)-Q(IQ0+1))*T7SA-(Q(IQ0+2)-Q(IQ0+1))*T7SC)/T7SE
      Q(1)    = Q(IQ0+1) - DQDU0*CS(IQ0+1)**2-.5*D2QDU20*CS(IQ0+1)**4
      DO I=2,IQ0-1
         Q(I) = Q(1) + DQDU0*CS(I)**2 + .5*D2QDU20*CS(I)**4
      ENDDO
       
C     COMPUTE P(1),DPDU0,D2PDU20, WHERE U=S^2
C     USING THREE POINTS PEQM(1),PEQ(2) AND PEQ(3)
      T7SA = (CSH(1)+CSH(2))**2 - CSH(1)**2
      T7SB = .5*((CSH(1)+CSH(2))**4-CSH(1)**4)
      T7SC = (CSH(1)+CSH(2)+CSH(3))**2 - CSH(1)**2
      T7SD = .5*((CSH(1)+CSH(2)+CSH(3))**4-CSH(1)**4)
      T7SE = T7SA*T7SD - T7SC*T7SB
      DPDU0   = ((PEQ(3)-PEQ(2))*T7SD-(PEQ(4)-PEQ(2))*T7SB)/T7SE
      D2PDU20 = ((PEQ(4)-PEQ(2))*T7SA-(PEQ(3)-PEQ(2))*T7SC)/T7SE
      PEQ(1)  = PEQ(2) - DPDU0*CSH(1)**2-.5*D2PDU20*CSH(1)**4
       
C     COMPUTE T7DUDR=(D U/D T7R), U=CS^2  
      DO I=1,NRP1
         T7DUDR(I) = T(I)/Q(I)*T7TMPS   
      ENDDO

C     COMPUTE T7DQDR=(D Q/D T7R) WITH SPLINE
      DO I=1,NR
         T7TMPA(I) = (Q(I+1)-Q(I))/(T7R(I+1)-T7R(I))
         T7TMPB(I) = (T7R(I+1)+T7R(I))*.5
      ENDDO
      CALL SPLINE1D(T7DQDR,T7R,NRP1,T7TMPA,T7TMPB,NR,WORKS(1,1))

C     COMPUTE T7D2QDR20=(D^2 Q/D T7R^2)(0)
      T7SA = D2QDU20*T7DUDR(1)
      T7SB = (TM(1)-T(1))/CSH(1)**2*4.*T7TMPS*DQDU0/Q(1)
      T7SC = -T7TMPS*T(1)*DQDU0**2/Q(1)**2
      T7D2QDR20 = (T7SA+T7SB+T7SC)*T7DUDR(1)

C     COMPUTE T7DPDR=(D P/D T7R) WITH SPLINE
      DO I=1,NR
         T7TMPA(I) = (PEQ(I+1)-PEQ(I))/(T7R(I+1)-T7R(I))
      ENDDO
      CALL SPLINE1D(T7DPDR,T7R,NRP1,T7TMPA,T7TMPB,NR,WORKS(1,1))

C     COMPUTE T7D2PDR20=(D^2 P/D T7R^2)(0)
      T7SA = D2PDU20*T7DUDR(1)
      T7SB = (TM(1)-T(1))/CSH(1)**2*4.*T7TMPS*DPDU0/Q(1)
      T7SC = -T7TMPS*T(1)*DQDU0*DPDU0/Q(1)**2
      T7D2PDR20 = (T7SA+T7SB+T7SC)*T7DUDR(1)

C     COMPUTE T7F FUNCTION
      DO I=1,NRP1
         T7F(I) = T7DUDR(I)*DPSIDS(NRP1)*.5/T7EPS0
      ENDDO

C     COMPUTE T7DGDR=(D G/D T7R)/T7F WITH SPLINE
      DO I=1,NR
         T7TMPA(I) = (T(I+1)-T(I))/(T7R(I+1)-T7R(I))
      ENDDO
      CALL SPLINE1D(T7DGDR,T7R,NRP1,T7TMPA,T7TMPB,NR,WORKS(1,1))
      DO I=1,NRP1
         T7DGDR(I) = T7DGDR(I)/T7F(I)
      ENDDO
      
C     RECOMPUTE T7DPDR=T7R*(D P/D T7R)/T7F^2
      DO I=1,NRP1
         T7DPDR(I) = T7DPDR(I)*T7R(I)/T7F(I)**2
      ENDDO

C     COMPUTE T7DFRDR=T7R^2*(D (T7F/T7R)/D R)/T7F WITH SPLINE
      DO I=1,NR
         T7TMPA(I) = (T7F(I+1)-T7F(I))/(T7R(I+1)-T7R(I))
      ENDDO
      CALL SPLINE1D(T7DFRDR,T7R,NRP1,T7TMPA,T7TMPB,NR,WORKS(1,1))
      DO I=1,NRP1
         T7DFRDR(I) = T7DFRDR(I)*T7R(I)/T7F(I) - 1.
      ENDDO

C     COMPUTE T7D2GDR2=T7R*(D (D G/D T7R)/T7F /D T7R) WITH SPLINE
      DO I=1,NR
         T7TMPA(I) = (T7DGDR(I+1)-T7DGDR(I))/(T7R(I+1)-T7R(I))
      ENDDO
      CALL SPLINE1D(T7D2GDR2,T7R,NRP1,T7TMPA,T7TMPB,NR,WORKS(1,1))
      DO I=1,NRP1
         T7D2GDR2(I) = T7D2GDR2(I)*T7R(I)
      ENDDO

C     COMPUTE RW2=1/|GRAD T7R|^2
      DO J=1,NCHI
      DO I=2,NRP1
         RW2(I,J) = (.5*RJA(I,J)*T7DUDR(I)/REQ(I,J)/CS(I))**2/G22L(I,J)
      ENDDO
      T7TMPS = (.5*RJAM(1,J)*T7DUDR(1)/REQM(1,J)/CSM(1))**2/G22LM(1,J)
      RW2(1,J) = (4*T7TMPS - RW2(2,J))/3.
      ENDDO

C     COMPUTE RW3=T7R*(GRAD T7R DOT GRAD CHI)/|GRAD T7R|^2
      DO J=1,NCHI
      DO I=2,NRP1
         RW3(I,J) = -.5*T7R(I)*T7DUDR(I)*G12L(I,J)/CS(I)/G22L(I,J)
      ENDDO
      T7TMPS = -.5*CSM(1)*G12LM(1,J)/G22LM(1,J)
      RW3(1,J) = (4*T7TMPS - RW3(2,J))/3.
      ENDDO
      
C     COMPUTE 2D EQUILIBRIUM QUANTITY: T7AK
      DO J=1,NCHI
      DO I=1,NRP1
         RW1(I,J) = REQ(I,J)**2
      ENDDO
      ENDDO                                                  
C
      call FFTDRIVER( RW1,  CW1,    FORWD, NRP1,  NRP1,   NPSTRT
     &                     ,MOUTP1, NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: CW1   in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      1,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW1,  CW1,    NRP1,  NRP1,   NPSTRT
     &                     ,MOUTP1, NCHI,  KUFFTP, 'CW1a')
      ENDIF
C
      DO I=1,NRP1
         T7AK(I,0) = DREAL(CW1(I,1))
      ENDDO
      DO M=1,MOUTPUT
      DO I=1,NRP1
         T7AK(I,M) = 2.*DREAL(CW1(I,M+1))
      ENDDO
      ENDDO

C     COMPUTE 2D EQUILIBRIUM QUANTITY: T7BK
      DO J=1,NCHI
      DO I=1,NRP1
         RW1(I,J) = RW2(I,J)/T7EPS0**2/REQ(I,J)**2
      ENDDO
      ENDDO                                                    
C
      call FFTDRIVER( RW1,  CW1,    FORWD, NRP1,  NRP1,   NPSTRT
     &                     ,MOUTP1, NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: CW1   in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      2,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW1,  CW1,    NRP1,  NRP1,   NPSTRT
     &                     ,MOUTP1, NCHI,  KUFFTP, 'CW1b')
      ENDIF
C
      DO I=1,NRP1
         T7BK(I,0) = DREAL(CW1(I,1))
      ENDDO
      DO M=1,MOUTPUT
      DO I=1,NRP1
         T7BK(I,M) = 2.*DREAL(CW1(I,M+1))
      ENDDO
      ENDDO

C     COMPUTE 2D EQUILIBRIUM QUANTITY: T7CK
      DO J=1,NCHI
      DO I=1,NRP1
         RW1(I,J) = RW2(I,J)/T7EPS0**2
      ENDDO
      ENDDO                                                    
C
      call FFTDRIVER( RW1,  CW1,    FORWD, NRP1,  NRP1,   NPSTRT
     &                     ,MOUTP1, NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: CW1   in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      3,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW1,  CW1,    NRP1,  NRP1,   NPSTRT
     &                     ,MOUTP1, NCHI,  KUFFTP, 'CW1c')
      ENDIF
C
      DO I=1,NRP1
         T7CK(I,0) = DREAL(CW1(I,1))
      ENDDO
      DO M=1,MOUTPUT
      DO I=1,NRP1
         T7CK(I,M) = 2.*DREAL(CW1(I,M+1))
      ENDDO
      ENDDO

C     COMPUTE 2D EQUILIBRIUM QUANTITY: T7DK
      DO J=1,NCHI
      DO I=1,NRP1
         RW1(I,J) = RW2(I,J)/T7EPS0**2*REQ(I,J)**2
      ENDDO
      ENDDO                                         
C
      call FFTDRIVER( RW1,  CW1,    FORWD, NRP1,  NRP1,   NPSTRT
     &                     ,MOUTP1, NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: CW1   in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      4,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW1,  CW1,    NRP1,  NRP1,   NPSTRT
     &                     ,MOUTP1, NCHI,  KUFFTP, 'CW1d')
      ENDIF
C
      DO I=1,NRP1
         T7DK(I,0) = DREAL(CW1(I,1))
      ENDDO
      DO M=1,MOUTPUT
      DO I=1,NRP1
         T7DK(I,M) = 2.*DREAL(CW1(I,M+1))
      ENDDO
      ENDDO

C     COMPUTE 2D EQUILIBRIUM QUANTITY: T7EK
      DO J=1,NCHI
      DO I=1,NRP1
         RW1(I,J) = RW2(I,J)/T7EPS0**2*REQ(I,J)**4
      ENDDO
      ENDDO                                                   
C
      call FFTDRIVER( RW1,  CW1,    FORWD, NRP1,  NRP1,   NPSTRT
     &                     ,MOUTP1, NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: CW1   in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      5,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW1,  CW1,    NRP1,  NRP1,   NPSTRT
     &                     ,MOUTP1, NCHI,  KUFFTP, 'CW1e')
      ENDIF
C
      DO I=1,NRP1
         T7EK(I,0) = DREAL(CW1(I,1))
      ENDDO
      DO M=1,MOUTPUT
      DO I=1,NRP1
         T7EK(I,M) = 2.*DREAL(CW1(I,M+1))
      ENDDO
      ENDDO

C     COMPUTE 2D EQUILIBRIUM QUANTITY: T7FK
      DO J=1,NCHI
      DO I=1,NRP1
         RW1(I,J) = RW3(I,J)
      ENDDO
      ENDDO                                              
C
      call FFTDRIVER( RW1,  CW1,    FORWD, NRP1,  NRP1,   NPSTRT
     &                     ,MOUTP1, NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: CW1   in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      6,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW1,  CW1,    NRP1,  NRP1,   NPSTRT
     &                     ,MOUTP1, NCHI,  KUFFTP, 'CW1f')
      ENDIF
C
      DO I=1,NRP1
         T7FK(I,0) = DIMAG(CW1(I,1))
      ENDDO
      DO M=1,MOUTPUT
      DO I=1,NRP1
         T7FK(I,M) =-2.*DIMAG(CW1(I,M+1))
      ENDDO
      ENDDO

C     COMPUTE 2D EQUILIBRIUM QUANTITY: T7GK
      DO J=1,NCHI
      DO I=1,NRP1
         RW1(I,J) = RW3(I,J)*REQ(I,J)**2
      ENDDO
      ENDDO                                                  
C
      call FFTDRIVER( RW1,  CW1,    FORWD, NRP1,  NRP1,   NPSTRT
     &                     ,MOUTP1, NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: CW1   in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      7,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW1,  CW1,    NRP1,  NRP1,   NPSTRT
     &                     ,MOUTP1, NCHI,  KUFFTP, 'CW1g')
      ENDIF
C
      DO I=1,NRP1
         T7GK(I,0) = DIMAG(CW1(I,1))
      ENDDO
      DO M=1,MOUTPUT
      DO I=1,NRP1
         T7GK(I,M) =-2.*DIMAG(CW1(I,M+1))
      ENDDO
      ENDDO

C     COMPUTE 2D EQUILIBRIUM QUANTITY: T7DADRK
      DO M=0,MOUTPUT
      DO I=1,NR
         T7TMPA(I) = (T7AK(I+1,M)-T7AK(I,M))/(T7R(I+1)-T7R(I))
      ENDDO
      CALL SPLINE1D(T7DC0DR,T7R,NRP1,T7TMPA,T7TMPB,NR,WORKS(1,1))
      DO I=1,NRP1
         T7DADRK(I,M) = T7DC0DR(I)*T7R(I)
      ENDDO
      ENDDO

C     COMPUTE T7D2GC0DR2=T7R*(D (D G/D T7R)*T7CK(0)/T7F /D T7R) WITH SPLINE
      DO I=2,NR
         T7TMPA(I) = (T7DGDR(I+1)*T7CK(I+1,0)*T7R(I+1)-
     &                T7DGDR(I)*T7CK(I,0)*T7R(I))
     &               /(T7R(I+1)-T7R(I))
      ENDDO
      CALL SPLINE1D(T7D2GC0DR2,T7R,NRP1,T7TMPA(2),T7TMPB(2),
     &              NR-1,WORKS(1,1))
      DO I=1,NRP1
         T7D2GC0DR2(I) = T7D2GC0DR2(I)-T7DGDR(I)*T7CK(I,0)
      ENDDO

C     COMPUTE T7D2PD0DR2=T7R*(D (D P/D T7R)*T7R*T7DK(0)/T7F^2 /D T7R) WITH SPLINE
      DO I=2,NR
         T7TMPA(I) = (T7DPDR(I+1)*T7DK(I+1,0)-T7DPDR(I)*T7DK(I,0))
     &               /(T7R(I+1)-T7R(I))
      ENDDO
      CALL SPLINE1D(T7D2PD0DR2,T7R,NRP1,T7TMPA(2),T7TMPB(2),
     &              NR-1,WORKS(1,1))
      DO I=1,NRP1
         T7D2PD0DR2(I) = T7D2PD0DR2(I)*T7R(I)
      ENDDO

C     COMPUTE T7DC0DR=T7R*(D T7CK(0)/D T7R)
      DO I=2,NR
         T7TMPA(I) = (T7CK(I+1,0)*T7R(I+1)-T7CK(I,0)*T7R(I))
     &               /(T7R(I+1)-T7R(I))
      ENDDO
      CALL SPLINE1D(T7DC0DR,T7R,NRP1,T7TMPA(2),T7TMPB(2),
     &              NR-1,WORKS(1,1))
      DO I=1,NRP1
         T7DC0DR(I) = T7DC0DR(I) - T7CK(I,0)
      ENDDO

C     OUTPUT ALL EQUILIBRIUM DATA INTO A FILE
      OPEN(CHOUTP,FILE='EQ4T7.OUT')
      REWIND(CHOUTP)
      WRITE(CHOUTP,1171) NRP1,MOUTPUT,T7EPS0,T7D2QDR20,T7D2PDR20,
     &                   DQDU0*T7DUDR(1),DPDU0*T7DUDR(1),0
      DO I=1,NRP1
         WRITE(CHOUTP,1172) CS(I),T7R(I),Q(I),T7DQDR(I),T7DGDR(I),
     &                      T7DPDR(I),T7DFRDR(I),T7D2GDR2(I)
      ENDDO
      DO I=1,NRP1
         WRITE(CHOUTP,1172) T7D2GC0DR2(I),T7D2PD0DR2(I),T7DC0DR(I),
     &                      T7F(I),T(I),PEQ(I),T7DUDR(I),0.0
      ENDDO
      DO M=0,MOUTPUT
      DO I=1,NRP1
         WRITE(CHOUTP,1172) T7AK(I,M),T7DADRK(I,M),T7BK(I,M),T7CK(I,M),
     &                      T7DK(I,M),T7EK(I,M),T7FK(I,M),T7GK(I,M)
      ENDDO
      ENDDO
      CLOSE(CHOUTP)
 1171 FORMAT(2(1X,I5),5(1X,E16.9),1(1X,I2))
 1172 FORMAT(E16.9,7(1X,E16.9))

      CS(1) = CS00

      DEALLOCATE(
     &   T7R,    T7TMPA,    T7DUDR,  T7DQDR,
     &   T7TMPB, T7DPDR,    T7F,     T7DGDR,
     &   T7DFRDR,T7D2GDR2,  T7D2GC0DR2, 
     &   T7D2PD0DR2,        T7DC0DR,
     &   RW1,      RW2,     RW3,
     &   WORKS,  CW1,
     &   T7AK,T7BK,T7CK,T7DK,T7EK,T7FK,T7GK,T7DADRK
     &        )

      RETURN
      END
