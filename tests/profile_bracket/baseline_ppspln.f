         SUBROUTINE BASELINE_PPSPLN(KN,PP,KPP,PS,RPP,D2RPP,PT,KD)
C        ############################################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
*  C2SP03  CUBIC SPLINE INTEPOLATION OF EXPERIMENTAL P-PRIME PROFILE  *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBAL.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSUR.inc'
C
         DIMENSION
     I   I1(NPT+2*NPISO),  IC(NPT+2*NPISO),
     R   D2RPP(NPISO),   PP(KN),   PS(NPISO),
     R   RPP(NPISO),     PT(KN)
C
C
C   BRACKET OUT [PS(I); PS(I+1)] INTERVAL SUCH THAT
C   PSIISO(I) <= PP(J) <= PSIISO(I+1), J=1,...,KN
C 
         IF (KPP.GT.NPISO) STOP 'KPP>NPISO'

         CALL RESETI(IC,KN,1)
         DO 1 JS = 1,KPP+1
           DO 1 JG=1,KN
             IF (IC(JG).EQ.0) GOTO 1
             ZS1 = 1. - PP(JG) / SPSIM
             IF (ZS1 .LT. 0.) ZS1 = 0.
             I1(JG) = JS-1
             IF (SQRT(ZS1).LE.PS(JS)) IC(JG) = 0
 1       CONTINUE
C
***********************************************************************
*                                                                     *
*  COMPUTE P-PRIME                                                    *
*                                                                     *
***********************************************************************
C
         DO 2 J2=1,KN
C
         IF (I1(J2) .LT. 1)   I1(J2) = 1
         IF (I1(J2) .GT. KPP) I1(J2) = KPP
C
         ZS1 = 1. - PP(J2) / SPSIM
C
         IF (ZS1 .LT. 0.) ZS1 = 0.
C
         ZS1 = SQRT(ZS1)
C
         ZH = PS(I1(J2)+1) - PS(I1(J2))
         ZA = (PS(I1(J2)+1) - ZS1) / ZH
         ZB = (ZS1 - PS(I1(J2))) / ZH
         ZC = (ZA + 1) * (ZA - 1) * ZH * (PS(I1(J2)+1) - ZS1) / 6.
         ZD = (ZB + 1) * (ZB - 1) * ZH * (ZS1 - PS(I1(J2))) / 6.
C 
CYQL2018
         IF (KD.EQ.0) THEN
         PT(J2) = ZA*RPP(I1(J2))   + ZB*RPP(I1(J2)+1) +
     +            ZC*D2RPP(I1(J2)) + ZD*D2RPP(I1(J2)+1)

         ELSE IF (KD.EQ.1) THEN
C        FIRST DERIVATIVE
         PT(J2) =-1./ZH*RPP(I1(J2)) + 1./ZH*RPP(I1(J2)+1) +
     +           (1./3.-ZA*ZA)*ZH/2.*D2RPP(I1(J2)) +
     +           (ZB*ZB-1./3.)*ZH/2.*D2RPP(I1(J2)+1)

         ENDIF

C
    2    CONTINUE
C         
         RETURN
         END
