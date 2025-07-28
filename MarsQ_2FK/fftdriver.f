      SUBROUTINE FFTDRIVER(RSP,FFT,ICASE,MDIM,M,MST,NF,NC
     &                            ,KUNIT,IERSUB,IERPLC,IERROR)
C
C     RSP     =   REAL SPACE FUNCTION    (INPUT FOR ICASE > 0, OUTPUT FOR ICASE < 0)
C     FFT     =   FAST FOURIER TRANSFORM (INPUT FOR ICASE < 0, OUTPUT FOR ICASE > 0)
C     ICASE   =   +1 FOR REAL TO FFT AND -1 FOR FFT TO REAL
C     MDIM    =   FIRST DIMENSION OF INPUT RSP OR FFT
C     M       =   NUMBER OF POINTS UP TO WHICH FFT OR ITS INVERSE IS DONE
C     MST     =   STARTING VALUE OF FIRST INDEX FOR INPUT RSP OR FFT
C     NF      =   ACTUAL NUMBER OF FFT VALUES IN FFT
C     NC      =   ACTUAL NUMBER OF ANGLES IN REAL SPACE
C     KUNIT   =   UNIT NUMBER FOR STANDARD OUTPUT
C     IERSUB  =   SUBROUTINE NAME FOR FLAGGED ERROR
C     IERPLC  =   ERROR LOCATION IDENTIFIER
C     IERROR  =   ERROR FLAG
C
C
      IMPLICIT NONE
C
      CHARACTER*16 IERSUB
      INTEGER      ICASE
      INTEGER      MDIM,     M,       NF,      NC,      MST
      INTEGER      KUNIT,    IERPLC,  IERROR,  INFLAG,  IEFLAG
      REAL*8       RSP(MDIM,NC)
      COMPLEX*16   FFT(MDIM,NF)
C
      INTEGER      INITL
C
      CHARACTER*16 SUBNAM,   ISFLAG
      INTEGER      IV,       IVP,      NV,       NPV,      NQV,
     &             NF1,      NF2,      NF2P,     NC2,      NCD,
     &             NDIMP,    MDIMP,    M0,       MSTP,     NK0
      INTEGER      NKFACTR
      REAL*8       ZERO
      REAL*8       NORMALIZE, RESLTIM
C
      INTEGER  ,   DIMENSION(:),   ALLOCATABLE::KFACTOR
      REAL*8,      DIMENSION(:),   ALLOCATABLE::PHASER
      REAL*8,      DIMENSION(:),   ALLOCATABLE::PHASEI
      REAL*8,      DIMENSION(:,:), ALLOCATABLE::RESULT
      REAL*8,      DIMENSION(:,:), ALLOCATABLE::STORAG
C
C
C
C 1.0 INITIALIZATION
C
C 1.1 INITIALIZE DATA
C
C 1.1.1 ERROR FLAGS
C
      SUBNAM  =  'FFTDRIVER'
      IERSUB  =  ''
      ISFLAG  =  'DEFAULT'
      IERROR  =  0
      IERPLC  =  0
      INFLAG  =  0
      IEFLAG  =  0
C
      ZERO    = 0.0D0
C
C 1.1.2 DIMENSIONS
C
      NF1     = NF - 1
      NF2     = 2*NF1
      NF2P    = NF2  + 1
      NC2     = NC/2 + 1
      NCD     = NC - NF + 2
C
      NK0     = NC
C
C 1.1.3 SET INITL TO FORCE THE PHASE ARRAYS TO BE CALCULATED
C
      INITL   =  0
C
C
C 1.2 CHECK INPUT STARTING VALUES AND DIMENSIONS
C
      NDIMP  =   NC
C
C
      IF(M .GT. MDIM) THEN
         IERSUB   =  SUBNAM
         IERPLC   = -1
         IERROR   =  MDIM - M
         RETURN
      ENDIF
C
      IF    (MST .GT. 0) THEN
         MSTP     =  MST
      ELSEIF(MST .LE. 0) THEN
         IERSUB   =  SUBNAM
         IERPLC   = -2
         IERROR   =  MST
         MSTP     = +1
      ENDIF
C
      MDIMP  =  M - MSTP + 1
      M0     =  MDIMP 
      IF(MDIMP .LT. 1) THEN
         IERSUB   =  SUBNAM
         IERPLC   = -3
         IERROR   =  MDIMP
         RETURN
      ENDIF
C
C
C 1.3 CHECK VALIDITY OF INPUT
C
C 1.3.1 ENSURE NC IS EVEN
C
      IF(2*(NC/2) .NE. NC) THEN
         IERSUB   =  SUBNAM
         IERPLC   = -4
         IERROR   =  NC
         RETURN
      ENDIF
C
C 1.3.2 ENSURE NF DOES NOT EXCEED NC/2 + 1
C
      IF(NF .GT. NC2) THEN
         IERSUB   =  SUBNAM
         IERPLC   = -5
         IERROR   =  NF
         RETURN
      ENDIF
C
C 1.3.3 ENSURE NF IS NON-NEGATIVE
C
      IF(NF .LE. 0) THEN
         IERSUB   =  SUBNAM
         IERPLC   = -6
         IERROR   = -NF
         RETURN
      ENDIF
C
C 1.3.3 SET A WARNING IF NF IS BELOW 2
C
      IF(NF .LE. 2) THEN
         IERSUB   =  SUBNAM
         IERPLC   = -7
         IERROR   = -NF
      ENDIF
C
C 1.3.1 SET A WARNING IF ICASE IS NOT CORRECTLY SET
C
      IF(IABS(ICASE) .NE. +1) THEN
         IERSUB   =  SUBNAM
         IERPLC   = -8
         IF(ICASE .EQ. 0) IERROR   = -1
         IF(ICASE .NE. 0) IERROR   =  ICASE
      ENDIF
C
C
C 1.4 ALLOCATE LOCAL ARRAYS
C
      ALLOCATE(KFACTOR(NC))
      ALLOCATE(PHASER(NDIMP))
      ALLOCATE(PHASEI(NDIMP))
      ALLOCATE(RESULT(MDIMP,NDIMP))
      ALLOCATE(STORAG(MDIMP,NDIMP))
C
C
C 1.5 INITIALIZE WORK ARRAY
C
      DO NV         = 1,NDIMP
      DO IV         = 1,MDIMP
      RESULT(IV,NV) =  ZERO
      STORAG(IV,NV) =  ZERO
      ENDDO
      ENDDO
C
C
C
C 2.0 INITIALIZE THE RESULT ARRAY WITH RSP OR FFT
C
C 2.1 STORE FROM 2 TO NF-1 AND FROM NF+1 TO 2(NF-1)
C
      IF(NF1 .GE. 2) THEN
        DO NV        = 2,NF1
        NPV          = NC - NV + 2
        NQV          = NF + NV - 1
        DO IV        = 1,MDIMP
        IVP          = IV + MSTP -1 
C
C     FORWARD FFT:
        IF    (ICASE .GT. 0) THEN
           RESULT(IV,NV)   =       RSP(IVP,NV)
           RESULT(IV,NQV)  =       RSP(IVP,NQV)
C
C     INVERSE FFT:
        ELSEIF(ICASE .LT. 0) THEN
           RESULT(IV,NV)   = +DREAL(FFT(IVP,NV))
           RESULT(IV,NPV)  = -DIMAG(FFT(IVP,NV))
        ENDIF
        ENDDO
        ENDDO
      ENDIF
C
C
C 2.2 STORE THE FIRST AND NF ELEMENT
C
      DO IV        = 1,MDIMP
      IVP          = IV + MSTP -1
      IF    (ICASE .GT. 0) THEN
            RESULT(IV, 1 )  =        RSP(IVP, 1 )
         IF(NF .GT. 1) THEN
            RESULT(IV, NF)  =        RSP(IVP, NF)
         ENDIF
C
      ELSEIF(ICASE .LT. 0) THEN
            RESULT(IV, 1 )  = +DREAL(FFT(IVP, 1 ))
         IF(NF .GT. 1) THEN
            RESULT(IV, NF)  = +DREAL(FFT(IVP, NF))
         ENDIF
         IF(NF .LT. NC2) THEN
            RESULT(IV, NCD) = -DIMAG(FFT(IVP, NF))
         ENDIF
      ENDIF
      ENDDO
C
C
C 2.3 FOR ICASE = +1 STORE THE REMAINING ELEMENTS FROM 2(NF-1) TO NC
C
      IF(ICASE .GT. 0) THEN
         IF(NF2P .LE. NC) THEN
            DO NV        = NF2P,NC
            DO IV        = 1,MDIMP
            IVP          = IV + MSTP -1
            RESULT(IV, NV)  =       RSP(IVP, NV)
            ENDDO
            ENDDO
         ENDIF
      ENDIF
C
C
C
C 3.0 CALCULATE THE FFT OR INVERSE FFT
C

      CALL FFTCALD(RESULT,STORAG,PHASER,PHASEI,M0,MDIMP,NC,NDIMP,
     &             KFACTOR,NKFACTR,NK0,ICASE,INITL,
     &             KUNIT,ISFLAG,INFLAG,IEFLAG)
C
      IF (IEFLAG .NE. 0) THEN
         IF    (INFLAG .EQ. 0) THEN
            WRITE(KUNIT,1000) INFLAG,IEFLAG
         ELSEIF(INFLAG .EQ. -3) THEN
            WRITE(KUNIT,2000) INFLAG,IEFLAG
         ELSEIF(INFLAG .EQ. -2) THEN
            WRITE(KUNIT,3000) INFLAG,IEFLAG
         ELSEIF(INFLAG .EQ. -1) THEN
            WRITE(KUNIT,4000) INFLAG,IEFLAG
         ELSEIF(INFLAG .GE. +1) THEN
            IF    (IEFLAG. EQ. 1) THEN
               WRITE(KUNIT,5000) M
            ELSEIF(IEFLAG .EQ. 2) THEN
               WRITE(KUNIT,5010) NC
            ELSEIF(IEFLAG .EQ. 3) THEN
               WRITE(KUNIT,5020) INITL
            ELSEIF(IEFLAG .EQ. 4) THEN
               WRITE(KUNIT,5030) INITL
            ELSEIF(IEFLAG .EQ. 5) THEN
               WRITE(KUNIT,5040) INITL
            ENDIF
         ENDIF
C
         IERSUB   = ISFLAG
         IERPLC   = INFLAG
         IERROR   = IEFLAG
         RETURN
      ENDIF
C
C
C
C 4.0 STORT THE RESULT AND CLEAN UP
C
C 4.1 RESET PARAMETERS FOR INITIALIZING PHASE ARRAYS
C
      INITL      = +1
C
C
C 4.2 STORE THE FINAL RESULT
C     THE RESULT IS STORED IN THE OUTPUT ARRAY FROM IV = 1,MDIMP AND SO
C     MAY BE MISALIGNED WITH THE INPUT ARRAY WHICH USES ONLY IV = MST,M
C
C 4.2.1 STORE THE NONZERO ELEMENTS
C
      IF(ICASE .GT. 0) NORMALIZE    = 1.0/DSQRT(DFLOAT(NC)) 
      IF(ICASE .LT. 0) NORMALIZE    =     DSQRT(DFLOAT(NC))
C
      IF(NF1 .GE. 2) THEN
        DO NV        = 2,NF1
        NPV          = NC - NV + 2
        NQV          = NF + NV - 1
        DO IV        = 1,MDIMP
        IVP          = IV + MSTP - 1
        IF    (ICASE .GT. 0) THEN
           FFT (IVP,NV) = NORMALIZE*CMPLX(RESULT(IV,NV),RESULT(IV,NPV))
        ELSEIF(ICASE .LT. 0) THEN
           RSP(IVP,NV)  = NORMALIZE*RESULT(IV,NV)
           RSP(IVP,NQV) = NORMALIZE*RESULT(IV,NQV)
        ENDIF
        ENDDO
        ENDDO
      ENDIF
C
      DO IV        = 1,MDIMP
      IVP          = IV + MSTP - 1
      IF(ICASE .GT. 0) THEN
            FFT(IVP, 1 )   = NORMALIZE*CMPLX(RESULT(IV, 1 ),   ZERO   )
         IF(NF .GT. 1  .AND.  NF .LE. NC2) THEN
            IF(NF .LT. NC2  .AND.  NCD .LE. NC) RESLTIM = RESULT(IV,NCD)
            IF(NF .EQ. NC2  .OR.   NCD .GT. NC) RESLTIM = ZERO
            FFT(IVP, NF)   = NORMALIZE*CMPLX(RESULT(IV, NF), RESLTIM  )
         ENDIF
C
      ELSEIF(ICASE .LT. 0) THEN
            RSP(IVP, 1 )   = NORMALIZE*RESULT(IV, 1 )
         IF(NF .GT. 1) THEN
            RSP(IVP, NF)   = NORMALIZE*RESULT(IV, NF)
         ENDIF
      ENDIF
      ENDDO
C
C
C 4.2.2 FILL UP TO THE DIMENSION
C
      IF(ICASE .LT. 0) THEN
         IF(NF2P .LE. NC) THEN
            DO NV        = NF2P,NC
            DO IV        = 1,MDIMP
            IVP          = IV + MSTP - 1
            RSP(IVP, NV)   = NORMALIZE*RESULT(IV, NV)
            ENDDO
            ENDDO
         ENDIF
      ENDIF
C
C
C 4.3 DEALLOCATE ARRAYS
C
      DEALLOCATE(KFACTOR)
      DEALLOCATE(PHASER)
      DEALLOCATE(PHASEI)
      DEALLOCATE(RESULT)
      DEALLOCATE(STORAG)
C
C
C
C 5.0 RETURN AND END
C
C     SUPPRESS ERROR MESSAGES
      IERROR = 0

      RETURN
C
 1000 FORMAT(1X,' ** ERROR IN FFTCALD:'
     &      ,1X,'INFLAG = ',I5,2X,'IEFLAG = ',I5)
 2000 FORMAT(1X,' ** ERROR IN PHASE:'
     &      ,1X,'INFLAG = ',I5,2X,'IEFLAG = ',I5)
 3000 FORMAT(1X,' ** ERROR IN FFTCAL0:'
     &      ,1X,'INFLAG = ',I5,2X,'IEFLAG = ',I5)
 4000 FORMAT(1X,' ** ERROR IN FFTCAL0:'
     &      ,1X,'INFLAG = ',I5,2X,'IEFLAG = ',I5)
 5000 FORMAT(1X,' ** M must be at least 1: M = ',I16)
 5010 FORMAT(1X,' ** N must be at least 1: N = ',I16)
 5020 FORMAT(1X,' ** ',I1,' is an invalid value of INITL')
 5030 FORMAT(1X,' ** INITL = ',I1,', but TRIG array never initialized')
 5040 FORMAT(1X,' ** INITL = ',I1,', but N and TRIG array incompatible')
      END
