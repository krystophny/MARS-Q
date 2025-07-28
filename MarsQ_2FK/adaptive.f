C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C GO THROUGH S-MESH ADAPTIVITY
C   1) FIND RATIONAL SURFACES
C   2) SHIFT S-MESH NEAR RATIONAL SURFACES
C   3) SPLINE INTERPOLATION OF ALL OUTRMAR QUANTITIES AT NEW MESH  
C
C YQ LIU, MARCH 20, 2012
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE GETADAPS
      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE ADAPTIVEM
      IMPLICIT NONE

C     GET NEW S-MESH
      CALL GETSMESH

C     SPLINE-INTERPOLATION OF Q ON NEW S-MESH
      IF (IADAPS.GT.0.OR.KADAPS.GT.0) CALL GETQPLSNEWMESH

C     COMPUTE ALL RATIONAL SURFACES
      CALL GETRATSURF

C     SHIFT S-MESH NEAR RATIONAL SURFACES
      IF (MRATSURF.GT.0) CALL SHIFTSRATSURF

C     SPLINE-INTERPOLATION OF ALL EQUILIBRIUM QUANTITIES ON NEW S-MESH
      IF (MRATSURF.GT.0.OR.IADAPS.GT.0.OR.KADAPS.GT.0) CALL GETEQNEWMESH

C     COMPUTE CHI(THETA) 
      IF (NCONVB1.EQ.1) THEN 
         CALL THETA2CHI
      ELSEIF (NCONVB1.EQ.2) THEN
         CALL B1EQAC2B1PEST
      ENDIF

      CS(1:NRP1)=CSE(1:NRP1)
      CSM(1:NR)=CSEM(1:NR)

      RETURN
      END

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C COMPUTE Q-PROFILE RIGHT AFTER READING EQUILIBRIUM DATA
C STORE Q-PROFILE IN <QPLS>
C YQ LIU, MARCH 20, 2012
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE GETQPLS
      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE GIJLM
      IMPLICIT NONE
      INCLUDE 'comioc.inc'
      INCLUDE 'comfft.inc'
C
      INTEGER NPSTRT,NUMFFT

      INTEGER I,J
      REAL*8,DIMENSION(:,:),ALLOCATABLE::RW1
      COMPLEX*16,DIMENSION(:,:),ALLOCATABLE::CTEMP

      IF (.NOT.ALLOCATED(RW1)) 
     &   ALLOCATE(RW1(NRP1,NCHI))
      IF (.NOT.ALLOCATED(QPLS)) ALLOCATE(QPLS(NRP1), QPLSM(NRP1))
      IF (.NOT.ALLOCATED(CTEMP)) ALLOCATE(CTEMP(NRP1,1))
C
      INCLUDE 'setfft.inc'
      SUBNAM    = 'GETQPLS'

      DO J=1,NCHI
         DO I=2,NRP1
            RW1(I,J)=RJA(I,J)*T(I)/REQ(I,J)**2/DPSIDS(I)
         ENDDO
         RW1(1,J)=RW1(2,J)
      ENDDO
C
      NPSTRT    =  1
      NUMFFT    =  1
      call FFTDRIVER( RW1,  CTEMP,   FORWD, NRP1,  NRP1,   NPSTRT
     &                     ,NUMFFT,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error for  CTEMP    in ',IERSUB
         call ABORTRUN
     &        (SUBNAM,      1,   MESSAGE
     &        ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW1,    CTEMP,     NRP1,    NRP1,    NPSTRT
     &                       ,NUMFFT,    NCHI,    KUFFTP, 'CTEMP1')
      ENDIF
C
      QPLS(1:NRP1) = DREAL(CTEMP(1:NRP1,1))

      DO J=1,NCHI
         DO I=1,NR
            RW1(I,J)=RJAM(I,J)*TM(I)/REQM(I,J)**2/DPSIDSM(I)
         ENDDO
      ENDDO
C
      call FFTDRIVER( RW1,  CTEMP,   FORWD, NRP1,  NR,     NPSTRT
     &                     ,NUMFFT,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error for  CTEMP    in ',IERSUB
         call ABORTRUN
     &        (SUBNAM,      2,   MESSAGE
     &        ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW1,    CTEMP,     NRP1,    NR,      NPSTRT
     &                       ,NUMFFT,    NCHI,    KUFFTP, 'CTEMP2')
      ENDIF
C
      QPLSM(1:NR) = DREAL(CTEMP(1:NR,1))

      QPLS(1) = (4.*QPLSM(1) - QPLS(2))/3.

      DEALLOCATE(RW1)
      DEALLOCATE(CTEMP)

      RETURN
      END

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C GET S-MESH
C YQ LIU, APRIL 17, 2012
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE GETSMESH
      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE ADAPTIVEM
      IMPLICIT NONE

      INTEGER I

      IF (.NOT.ALLOCATED(CSEN))  ALLOCATE(CSEN(NRP1),CSENM(NRP1))
      IF (.NOT.ALLOCATED(QPLSO)) ALLOCATE(QPLSO(NRP1),QPLSOM(NRP1))

      QPLSO(1:NRP1)  = QPLS(1:NRP1)
      QPLSOM(1:NRP1) = QPLSM(1:NRP1)

      IF (IADAPS.EQ.0.AND.KADAPS.EQ.0) THEN
         CSEN(1:NRP1)  = CSE(1:NRP1)
      ENDIF

      IF (IADAPS.EQ.0.AND.KADAPS.EQ.1) THEN
         OPEN(99,FILE='PROFEQ.OUT')
         DO I=1,NRP1
            READ(99,120) CSEN(I)
         ENDDO
 120     FORMAT(16(E16.9,1X))
         CLOSE(99)
      ENDIF
         
      DO I=1,NR
         CSENM(I) = 0.5*(CSEN(I)+CSEN(I+1))
      ENDDO
      CSENM(NRP1) = CSEN(NRP1)

      RETURN
      END

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C SPLINE INTERPOLATION OF Q-PROFILE AT NEW MESH  
C YQ LIU, MARCH 20, 2012
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE GETQPLSNEWMESH
      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE GIJLM
      USE ADAPTIVEM
      IMPLICIT NONE

      REAL*8,DIMENSION(:),ALLOCATABLE::PSIISO,RW1

      ALLOCATE(RW1(NRP1),PSIISO(NRP1))

      CALL SPLINE1D(RW1,CSEN,NRP1,QPLSO,CSE,NRP1,PSIISO)
      QPLS(1:NRP1) = RW1(1:NRP1)         
      CALL SPLINE1D(RW1,CSENM,NR,QPLSOM,CSEM,NR,PSIISO)
      QPLSM(1:NR) = RW1(1:NR)         

      DEALLOCATE(RW1,PSIISO)

      RETURN
      END
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C FIND ALL RATIONAL SURFACES
C YQ LIU, MARCH 20, 2012
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE GETRATSURF
      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE GIJLM
      USE ADAPTIVEM
      USE MPIENV
      IMPLICIT NONE
      INTEGER I,J,K,MMIN,MMAX
      REAL*8  QMIN,QMAX

C     FIND TOTAL NUMBER OF RATIONAL SURFACES
C     FOR GENERIC (EVEN NON-MONOTONIC) Q-PROFILE
      QMIN=QPLS(1)
      QMAX=QPLS(NRP1)
      DO I=1,NRP1
         IF (QPLS(I).LT.QMIN) QMIN=QPLS(I)
         IF (QPLS(I).GT.QMAX) QMAX=QPLS(I)
      ENDDO
     
      IF (RNTOR.LT.0) THEN
         MMIN     = MAX(INT(-RNTOR*QMIN)+1,M1)
         MMAX     = MIN(INT(-RNTOR*QMAX),M2)
      ELSE
         MMIN     = MAX(INT(-RNTOR*QMAX),M1)
         MMAX     = MIN(INT(-RNTOR*QMIN)-1,M2)
      ENDIF

      NRATSURF = 0
      DO I=1,NR
         DO J=MMIN,MMAX
            QMIN=(-QPLS(I)*RNTOR-J)*(-QPLS(I+1)*RNTOR-J)
            IF (QMIN.LE.0.) NRATSURF=NRATSURF+1
         ENDDO
      ENDDO

C     ALLOCATE ARRAYS
      IF (NRATSURF.GT.0) THEN
      IF (.NOT.ALLOCATED(IRATSURF)) ALLOCATE(IRATSURF(NRATSURF))
      IF (.NOT.ALLOCATED(QRATSURF)) ALLOCATE(QRATSURF(NRATSURF))
      ENDIF

C     IDENTIFY ALL RATIONAL SURFACES
      K = 0
      DO I=1,NR
         DO J=MMIN,MMAX
            QMIN=(-QPLS(I)*RNTOR-J)*(-QPLS(I+1)*RNTOR-J)
            IF (QMIN.LE.0.) THEN
               K = K+1
               IRATSURF(K) = I
               QRATSURF(K) = DFLOAT(J)/(-RNTOR)
            ENDIF
         ENDDO
      ENDDO
      
      RETURN
      END

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C SHIFT RADIAL MESH NEAR RATIONAL SURFACES ACCORDING TO <MRATSURF>
C    MRATSURF = 0: DO NOT SHIFT RADIAL MESH
C    MRATSURF = 1: INTEGER MESH TO ALIGN EXACTLY WITH RATIONAL SURFACE
C    MRATSURF = 2: HALF-INTEGER MESH TO ALIGN EXACTLY WITH RATIONAL SURFACE
C YQ LIU, MARCH 20, 2012
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE SHIFTSRATSURF
      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE GIJLM
      USE ADAPTIVEM
      IMPLICIT NONE

      INTEGER I,K
      REAL*8 H1, H2, X1, X2, X3, Y1, Y2, Y3
      REAL*8,DIMENSION(:),ALLOCATABLE::CSEO,CSEOM

      IF (.NOT.ALLOCATED(CSEO)) ALLOCATE(CSEO(NRP1),CSEOM(NRP1))

      CSEO(1:NRP1)  = CSEN(1:NRP1)
      CSEOM(1:NRP1) = CSENM(1:NRP1)
      
      DO K=1,NRATSURF
         I = IRATSURF(K)
         X1 = CSEO(I)
         X2 = CSEO(I+1)
         Y1 = QPLS(I)
         Y2 = QPLS(I+1)
         Y3 = QRATSURF(K)
         
         X3 = X1 + (Y3-Y1)/(Y2-Y1)*(X2-X1)

         IF (MRATSURF.EQ.1) THEN
            H1 = X3-X1
            H2 = X2-X3
            IF (H1.LE.H2.AND.I.GT.1) THEN
               CSEN(I) = X3
               CSENM(I) = (X3+X2)*0.5
               CSENM(I-1) = (CSEO(I-1)+X3)*0.5
            ENDIF
            IF (H1.GT.H2.AND.I.LT.NR) THEN
               CSEN(I+1)  = X3
               CSENM(I)   = (X1+X3)*0.5
               CSENM(I+1) = (X3+CSEO(I+2))*0.5
               IRATSURF(K)= I+1
            ENDIF
         ELSEIF (MRATSURF.EQ.2) THEN
            H2 = ((X1+X2)*0.5-X3)*0.5
            X1 = X1-H2
            X2 = X2-H2
            IF (I.GT.1.AND.I.LT.NR.AND.X1.GT.CSEO(I-1).AND.
     &          X2.LT.CSEO(I+2)) THEN
               CSEN(I)    = X1
               CSEN(I+1)  = X2
               CSENM(I-1) = (X1+CSEO(I-1))*0.5
               CSENM(I)   = X3
               CSENM(I+1) = (X2+CSEO(I+2))*0.5
            ENDIF
         ENDIF
      ENDDO

      DEALLOCATE(CSEO,CSEOM)

      RETURN
      END


C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C SPLINE INTERPOLATION OF ALL OUTRMAR QUANTITIES AT NEW MESH  
C YQ LIU, MARCH 20, 2012
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE GETEQNEWMESH
      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE GIJLM
      USE ADAPTIVEM
      IMPLICIT NONE

      INTEGER I,L,J
      REAL*8,DIMENSION(:),ALLOCATABLE::PSIISO,RW1

      ALLOCATE(RW1(NRP1),PSIISO(NRP1))

      CALL SPLINE1D(RW1,CSEN,NRP1,QPLSO,CSE,NRP1,PSIISO)
      QPLS(1:NRP1) = RW1(1:NRP1)         
      CALL SPLINE1D(RW1,CSENM,NR,QPLSOM,CSEM,NR,PSIISO)
      QPLSM(1:NR) = RW1(1:NR)         

      CALL SPLINE1D(RW1,CSEN,NRP1,PEQ,CSE,NRP1,PSIISO)
      PEQ(1:NRP1) = RW1(1:NRP1)         
      CALL SPLINE1D(RW1,CSENM,NR,PEQM,CSEM,NR,PSIISO)
      PEQM(1:NR) = RW1(1:NR)         

      CALL SPLINE1D(RW1,CSEN,NRP1,T,CSE,NRP1,PSIISO)
      T(1:NRP1) = RW1(1:NRP1)         
      CALL SPLINE1D(RW1,CSENM,NR,TM,CSEM,NR,PSIISO)
      TM(1:NR) = RW1(1:NR)         

      CALL SPLINE1D(RW1,CSEN,NRP1,PPEQ,CSE,NRP1,PSIISO)
      PPEQ(1:NRP1) = RW1(1:NRP1)         
      CALL SPLINE1D(RW1,CSENM,NR,PPEQM,CSEM,NR,PSIISO)
      PPEQM(1:NR) = RW1(1:NR)         

      CALL SPLINE1D(RW1,CSEN,NRP1,DPSIDS,CSE,NRP1,PSIISO)
      DPSIDS(1:NRP1) = RW1(1:NRP1)         
      CALL SPLINE1D(RW1,CSENM,NR,DPSIDSM,CSEM,NR,PSIISO)
      DPSIDSM(1:NR) = RW1(1:NR)         

      CALL SPLINE1D(RW1,CSEN,NRP1,TP,CSE,NRP1,PSIISO)
      TP(1:NRP1) = RW1(1:NRP1)         
      CALL SPLINE1D(RW1,CSENM,NR,TPM,CSEM,NR,PSIISO)
      TPM(1:NR) = RW1(1:NR)         

      DO L=1,NCHI
         CALL SPLINE1D(RW1,CSEN,NRP1,REQ(1,L),CSE,NRP1,PSIISO)
         REQ(1:NRP1,L) = RW1(1:NRP1)         
         CALL SPLINE1D(RW1,CSENM,NR,REQM(1,L),CSEM,NR,PSIISO)
         REQM(1:NR,L) = RW1(1:NR)         

         CALL SPLINE1D(RW1,CSEN,NRP1,ZEQ(1,L),CSE,NRP1,PSIISO)
         ZEQ(1:NRP1,L) = RW1(1:NRP1)         
         CALL SPLINE1D(RW1,CSENM,NR,ZEQM(1,L),CSEM,NR,PSIISO)
         ZEQM(1:NR,L) = RW1(1:NR)         

         CALL SPLINE1D(RW1,CSEN,NRP1,G11L(1,L),CSE,NRP1,PSIISO)
         G11L(1:NRP1,L) = RW1(1:NRP1)         
         CALL SPLINE1D(RW1,CSENM,NR,G11LM(1,L),CSEM,NR,PSIISO)
         G11LM(1:NR,L) = RW1(1:NR)         

         CALL SPLINE1D(RW1,CSEN,NRP1,G22L(1,L),CSE,NRP1,PSIISO)
         G22L(1:NRP1,L) = RW1(1:NRP1)         
         CALL SPLINE1D(RW1,CSENM,NR,G22LM(1,L),CSEM,NR,PSIISO)
         G22LM(1:NR,L) = RW1(1:NR)         

         CALL SPLINE1D(RW1,CSEN,NRP1,G33L(1,L),CSE,NRP1,PSIISO)
         G33L(1:NRP1,L) = RW1(1:NRP1)         
         CALL SPLINE1D(RW1,CSENM,NR,G33LM(1,L),CSEM,NR,PSIISO)
         G33LM(1:NR,L) = RW1(1:NR)         
 
         CALL SPLINE1D(RW1,CSEN,NRP1,G12L(1,L),CSE,NRP1,PSIISO)
         G12L(1:NRP1,L) = RW1(1:NRP1)         
         CALL SPLINE1D(RW1,CSENM,NR,G12LM(1,L),CSEM,NR,PSIISO)
         G12LM(1:NR,L) = RW1(1:NR)         

         CALL SPLINE1D(RW1,CSEN,NRP1,RDCDZ(1,L),CSE,NRP1,PSIISO)
         RDCDZ(1:NRP1,L) = RW1(1:NRP1)         
         CALL SPLINE1D(RW1,CSENM,NR,RDCDZM(1,L),CSEM,NR,PSIISO)
         RDCDZM(1:NR,L) = RW1(1:NR)         

         CALL SPLINE1D(RW1,CSEN,NRP1,RDSDZ(1,L),CSE,NRP1,PSIISO)
         RDSDZ(1:NRP1,L) = RW1(1:NRP1)         
         CALL SPLINE1D(RW1,CSENM,NR,RDSDZM(1,L),CSEM,NR,PSIISO)
         RDSDZM(1:NR,L) = RW1(1:NR)         

         CALL SPLINE1D(RW1,CSEN,NRP1,RBZ(1,L),CSE,NRP1,PSIISO)
         RBZ(1:NRP1,L) = RW1(1:NRP1)         
         CALL SPLINE1D(RW1,CSENM,NR,RBZM(1,L),CSEM,NR,PSIISO)
         RBZM(1:NR,L) = RW1(1:NR)         
      ENDDO

      CSE(1:NRP1)  = CSEN(1:NRP1)
      CSEM(1:NRP1) = CSENM(1:NRP1)

C     RE-COMPUTE JACOBIAN
      DO I= 1,NRP1
         DO L=1,NCHI
            RJA(I,L) = SQRT(ABS(G33L(I,L)*(G11L(I,L)*G22L(I,L)-
     &                      G12L(I,L)**2)))
         ENDDO
      ENDDO
      DO I= 1,NR
         DO L=1,NCHI
            RJAM(I,L) = SQRT(ABS(G33LM(I,L)*(G11LM(I,L)*G22LM(I,L)-
     &                       G12LM(I,L)**2)))
         ENDDO
      ENDDO

C     PERFORM OTHER ADJUSTMENTS 
      CSE(1) = 0.
      PEQ(1) = (4.*PEQM(1)-PEQ(2))/3.
      PPEQ(1) = 0.
      T(1)   = (4.*TM(1)  -T(2))  /3.
      TP (1) = 0.
      DPSIDS(1) = 0.
      REQ(1,1) = (4.*REQM(1,1)-REQ(2,1))/3.
      ZEQ(1,1) = (4.*ZEQM(1,1)-ZEQ(2,1))/3.
      TM(NRP1) = T(NRP1)
C
      DO 40 J = 1,NCHI
        REQ(1,J) = REQ(1,1)
        ZEQ(1,J) = ZEQ(1,1)
        RJA(1,J) = 0.
C       G11L(1,J) = 0.
        G22L(1,J) = 0.
C       G33L(1,J) = 0.
        G12L(1,J) = 0.
 40   CONTINUE

      DEALLOCATE(RW1,PSIISO)

      RETURN
      END

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C GET NEW S-MESH BASED ON RADIAL PROFILE OF COMPUTED TOROIDAL JXB TORQUE
C YQ LIU, APRIL 17, 2012
C KOPTION = 1: OLD METHOD TO GENERATE NEW MESH
C           2: USE LSODE TO GENERATE NEW MESH
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE GETSNEWMESH
      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE TORQUEM
      USE ADAPTIVEM
      USE TOOLBOX
      IMPLICIT NONE

      INTEGER I,J,NRP2,KCHECK
      REAL    TMP,TMP1,TMP2,TMP3
      REAL*8,DIMENSION(:),ALLOCATABLE::ADS,ADDL,ADD2,ADLF,ADLU
      REAL*8,DIMENSION(:,:),ALLOCATABLE::ADF

      KCHECK  = 1

      IF (ADAPSOPT.EQ.0) THEN
         ALLOCATE(ADS(NR),ADF(1,NR))

         ADS = CSM(1:NR)
         ADF(1,1:NR) = TORQUEJXB(1:NR)
C        ADF(1,1:NR) = 0.25-(CSM(1:NR)-0.5)**2

         CALL ADAPTIVE_GRID_LSODE(ADS,ADF)

         CSEN(1) = 0.
         DO I=1,NR-1
            CSEN(I+1)= 2.*ADS(I) - CSEN(I)
         ENDDO
         CSEN(NRP1) = 1.

         IF (KCHECK.EQ.1) THEN
            WRITE(*,*) 'ADAPTIVE_GRID_LSODE:'
            DO I=1,NR
               WRITE(*,110) CSM(I),ADF(1,I),ADS(I)
            ENDDO
         ENDIF 
 110     FORMAT(3(E16.9,1X))

         DEALLOCATE(ADS, ADF)
      
         RETURN
      ENDIF

      IF (ADAPSOPT.NE.1) STOP 'NO OPTION FOR ADAPTIVE PACKING'
      NRP2    = NRP1 + 1

      ALLOCATE(ADS(NRP2),  ADF(1,NRP2),  ADDL(NRP2),
     &         ADD2(NRP2), ADLF(NRP2), ADLU(NRP1))

C     TEMPARY S-MESH
      ADS(2:NRP1) = CSM(1:NR)
      ADS(1)      = 0.0
      ADS(NRP2)   = 1.0

C     INDICATOR FUNCTION (NORMALISED TORQUE DENSITY)
      ADF(1,2:NRP1) = TORQUEJXB(1:NR)
      ADF(1,1)      = ADF(1,2)
      ADF(1,NRP2)   = ADF(1,NRP1)
   
      TMP = 0.0
      DO I=1,NRP2
         IF (ABS(ADF(1,I)).GT.TMP) TMP = ABS(ADF(1,I))
      ENDDO
      IF (TMP.EQ.0.) TMP=1.0
      ADF = ADF/TMP*ADAPSDF

C     COMPUTE LENGTH OF EACH SEGMENT OF CURVE ADF(ADS)
      ADDL(1) = 0.0
      DO I=2,NRP2
         ADDL(I)=SQRT((ADS(I)-ADS(I-1))**2+(ADF(1,I)-ADF(1,I-1))**2)*
     &           Q(I-1)
      ENDDO

      TMP = 0.0
      DO I=1,NRP2
         IF (ABS(ADDL(I)).GT.TMP) TMP = ABS(ADDL(I))
      ENDDO
      ADDL = ADDL/TMP
      
C     COMPUTE INVERSE CURVATURE OF CURVE ADF(ADS) (NORMALISED SECOND DERIVATIVE)
      DO I=2,NRP1
         TMP1 = (ADF(1,I)-ADF(1,I-1))/(ADS(I)-ADS(I-1))              
         TMP2 = (ADF(1,I+1)-ADF(1,I))/(ADS(I+1)-ADS(I))
         TMP3 = 0.5*(ADS(I+1)-ADS(I-1))
         ADD2(I) = (TMP2-TMP1)/TMP3*Q(I-1)
      ENDDO
      ADD2(1)    = 0.0
      ADD2(NRP2) = ADD2(NRP1)
      ADD2       = ABS(ADD2)

      TMP = 0.0
      DO I=1,NRP2
         IF (ABS(ADD2(I)).GT.TMP) TMP = ABS(ADD2(I))
      ENDDO
      IF (TMP.EQ.0.) TMP=1.0
      ADD2 = ADD2/TMP

C     COMPUTE CRITERION FUNCTION FOR NEW MESH
      ADLF(1) = 0.0
      DO I=2,NRP2
         ADLF(I) = ADLF(I-1) + (ADDL(I)+ADD2(I)*ADAPSD2)**ADAPSLF
      ENDDO

      TMP  = ADLF(NRP2)
      ADLF = ADLF/TMP

C     SMOOTH ADLF
      TMP=1.0
      DO I=1,NRP1
         IF ((ADS(I+1)-ADS(I)).LT.TMP) TMP = ADS(I+1) - ADS(I)
      ENDDO

      DO J=1,20
         ADDL = ADLF
         DO I=2,NRP1
C           TMP1 = (ADDL(I)-ADDL(I-1))/(ADS(I)-ADS(I-1))              
C           TMP2 = (ADDL(I+1)-ADDL(I))/(ADS(I+1)-ADS(I))
C           TMP3 = 0.5*(ADS(I+1)-ADS(I-1))
C           ADLF(I) = ADDL(I) + (TMP2-TMP1)/TMP3*TMP**2*ADAPSNU
            ADLF(I) = ADDL(I-1)/3 + ADDL(I)/3 + ADDL(I+1)/3
         ENDDO
      ENDDO

C     DEFINE UNIFORM MESH FOR CRITERION FUNCTION
      TMP = 1.0/DFLOAT(NR)
      ADLU(1) = 0.0
      DO I=2,NR
         ADLU(I) = ADLU(I-1) + TMP
      ENDDO
      ADLU(NRP1) = 1.0

C     USE LINEAR ITERPOLATION TO FIND NEW S-MESH
      CSEN(1)    = 0.0
      CSEN(NRP1) = 1.0
      DO I=2,NR
         TMP  = ADLU(I)
         DO J=1,NRP1
            TMP1 = ADLF(J)
            TMP2 = ADLF(J+1)
            IF ((TMP-TMP1)*(TMP2-TMP).GE.0.0) 
     &      CSEN(I) = ADS(J) + (ADS(J+1)-ADS(J))*(TMP-TMP1)/(TMP2-TMP1)
         ENDDO
      ENDDO

C     CHECK MESH ADAPTIVITY 
      IF (KCHECK.EQ.1) THEN
         WRITE(*,*) 'ADAPTIVE NEW MESH:'
         DO I=1,NRP1
            WRITE(*,120) ADS(I),ADF(1,I),ADDL(I),ADD2(I),ADLF(I),ADLU(I)
     &                   ,CSEN(I)
         ENDDO
         I = NRP2
         WRITE(*,120) ADS(I),ADF(1,I),ADDL(I),ADD2(I),ADLF(I),0.0,0.0 
 120     FORMAT(16(E12.5,1X))
      ENDIF

      DEALLOCATE(ADS,ADF,ADDL,ADD2,ADLF,ADLU)

      RETURN
      END

