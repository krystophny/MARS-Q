      SUBROUTINE FFTCALD(RESULT,STORAG,PHASER,PHASEI,M0,MDIM,N0,NDIM
     &                  ,KFACTR,NKFACT,NK0,ICASE,INITL,KUNIT
     &                  ,IESUBR,INFLAG,IEFLAG)
C
      IMPLICIT NONE
C
      INTEGER      KUNIT
      INTEGER      NKFACT, NK0,    NKF
      INTEGER      KFACTR(NK0)
C
      CHARACTER*16 IESUBR
      INTEGER      M0,     N0,     MDIM,   NDIM
      INTEGER      ICASE
      INTEGER      IEFLG0, IEFLG1, INFLG0, INFLG1
C
      CHARACTER*16 SUBNAM, IESUB0, IESUB1
      INTEGER      INFLAG, IEFLAG
      INTEGER      INITL
      REAL*8       RESULT(MDIM,NDIM),STORAG(MDIM,NDIM)
      REAL*8       PHASER(NDIM),     PHASEI(NDIM)
C
C
C
C 1.0 INITIALIZATION
C
C 1.1 INITIALIZE ERROR FLAGS
C
      SUBNAM   = 'FFTCALD'
      IESUBR   = ''
      INFLAG   = 0
      IEFLAG   = 0
      IEFLG0   = 0
      IEFLG1   = 0
C
C
C 1.2 MISCELLANEOUS DATA
C
      NKF      = NK0
C
C
C 1.3 CHECK FOR INVALID INPUT
C
      IF(M0 .LT. 1) THEN
         IESUBR  = SUBNAM
         INFLAG  = +1
         IEFLAG  = +1
         RETURN
      ENDIF
C
      IF(N0 .LT. 1) THEN
         IESUBR  = SUBNAM
         INFLAG  = +2
         IEFLAG  = +1
         RETURN
      ENDIF
C
C
C 2.0 CALCULATE PHASE FACTORS FOR EACH CASE DEPENDING ON FACTORS OF N0
C
      CALL PHASE(PHASER,PHASEI,KFACTR,NKFACT,NKF,N0,NDIM,INITL
     &                               ,KUNIT,IESUB0,INFLG0,IEFLG0)
C
      IF(IEFLG0 .NE. 0) THEN
         IESUBR   =  IESUB0
         INFLAG   = +INFLG0
         IEFLAG   = +IEFLG0
         RETURN
      ENDIF   
C
C
C
C 3.0 CALCULATE FFT
C
      CALL FFTCAL0(RESULT,STORAG,PHASER,PHASEI,KFACTR,NKFACT
     &                          ,NKF,M0,MDIM,N0,NDIM,ICASE
     &                          ,KUNIT,IESUB1,INFLG1,IEFLG1)
C
      IF(IEFLG1 .NE. 0) THEN
         IESUBR   =  IESUB1
         INFLAG   = +INFLG1
         IEFLAG   = +IEFLG1
         RETURN
      ENDIF
C
C
C
C 4.0 RETURN AND END
C
      RETURN
      END
      SUBROUTINE FFTCAL0(RESULT,STORAG,PHASER,PHASEI,KFACTR,NKFACT
     &                             ,NKF,M0,MDIM,N0,NDIM,ICASE
     &                             ,KUNIT,IESUBR,IERPLC,IERROR)

C
C
      IMPLICIT NONE
C
      CHARACTER*16 SUBNAM, IESUB
      INTEGER      IV,     K,      KP,     KV,     KT
      INTEGER      NPAR,   KPAR
      INTEGER      ND1,    ND2,    ND3,    NST
C
      CHARACTER*16 IESUBR
      INTEGER      KUNIT
      INTEGER      ICASE
      INTEGER      M0,     N0,     MDIM,   NDIM
      INTEGER      NKF,    NKFACT
      INTEGER      KFACTR(NKF)
      INTEGER      IERPLC, IERROR
      INTEGER      MP,     NP,     NPD2,   MSIGN
C
      REAL*8       FACTOR
      REAL*8       RESULT(MDIM,NDIM), STORAG(MDIM,NDIM)
      REAL*8       PHASER(NDIM),      PHASEI(NDIM)
C
C
C
C 1.0 INITIALIZATION
C
C 1.1 INITIALIZE PARAMETERS
C
      MP      = M0
      NP      = N0
      NPD2    = NP/2 + 1
C
      SUBNAM  =  'FFTCAL0'
      IESUB   =  ''
      IESUBR  =  SUBNAM
      IERPLC  =  0
      IERROR  =  0
C
      MSIGN   = +1
      IF    (ICASE .GT. 0) THEN
         KPAR    = +1 
         NPAR    =  N0
      ELSEIF(ICASE .LT. 0) THEN
         KPAR    =  N0
         NPAR    = +1 
      ENDIF
C
C
C 1.2 CHECK FOR VALID INPUT
C
      IF(NKFACT .LT. 1) THEN
         IESUBR  =  SUBNAM
         IERPLC  =  1
         IERROR  =  NKFACT
         RETURN
      ENDIF
C
C
C 1.2 SET UP FOR INVERSE FFT
C
      IF    (ICASE .LT. 0) THEN
         DO IV        = 1, MP
         RESULT(IV,1) = 0.5*RESULT(IV,1)
         ENDDO
C
         IF(2*(NP/2) .EQ. NP) THEN
            DO IV           = 1, MP
            RESULT(IV,NPD2) = 0.5*RESULT(IV,NPD2)
            ENDDO
         END IF
      END IF
C
C
C
C 2.0 LOOP OVER KFACTR VALUES
C
      DO K       = 1,NKFACT
         KP         = K
         KV         = NKFACT - K + 1
         IF(ICASE .GT. 0) KT         = KFACTR(KV)
         IF(ICASE .LT. 0) KT         = KFACTR(KP)
         IF(KT .EQ. 0) THEN
            IESUBR  =  SUBNAM
            IERPLC  = +2
            IERROR  = +K
            RETURN
         ENDIF
C
         IF(ICASE .GT. 0) NPAR       =  NPAR/KT
         IF(ICASE .LT. 0) KPAR       =  KPAR/KT
         ND1        =  NPAR*M0
         ND2        =  KT
         ND3        =  KPAR
         NST        =  (NPAR-1)*KT*KPAR + 1
C
         IF    (ICASE .GT. 0) THEN
            IF    (KT .EQ. 2) THEN
               IESUBR   = 'FFTCAL2F'
               CALL FFTCAL2F(RESULT,STORAG,PHASER,PHASEI,MDIM,M0,NDIM,N0
     &                                           ,ND1,ND2,ND3,NST,MSIGN
     &                                           ,IESUB,IERPLC,IERROR)
            ELSEIF(KT .EQ. 3) THEN
               IESUBR   = 'FFTCAL3F'
               CALL FFTCAL3F(RESULT,STORAG,PHASER,PHASEI,MDIM,M0,NDIM,N0
     &                                           ,ND1,ND2,ND3,NST,MSIGN
     &                                           ,IESUB,IERPLC,IERROR)
            ELSEIF(KT .EQ. 4) THEN
               IESUBR   = 'FFTCAL4F'
               CALL FFTCAL4F(RESULT,STORAG,PHASER,PHASEI,MDIM,M0,NDIM,N0
     &                                           ,ND1,ND2,ND3,NST,MSIGN
     &                                           ,IESUB,IERPLC,IERROR)
            ELSEIF(KT .EQ. 5) THEN
               IESUBR   = 'FFTCAL5F'
               CALL FFTCAL5F(RESULT,STORAG,PHASER,PHASEI,MDIM,M0,NDIM,N0
     &                                           ,ND1,ND2,ND3,NST,MSIGN
     &                                           ,IESUB,IERPLC,IERROR)
            ELSEIF(KT .EQ. 6) THEN
               IESUBR   = 'FFTCAL6F'
               CALL FFTCAL6F(RESULT,STORAG,PHASER,PHASEI,MDIM,M0,NDIM,N0
     &                                           ,ND1,ND2,ND3,NST,MSIGN
     &                                           ,IESUB,IERPLC,IERROR)
            ELSEIF(KT .GT. 6) THEN
               IESUBR   = 'FFTCAL7F'
               CALL FFTCAL7F(RESULT,STORAG,PHASER,PHASEI,MDIM,M0,NDIM,N0
     &                                           ,ND1,ND2,ND3,NST,MSIGN
     &                                           ,IESUB,IERPLC,IERROR)
            ENDIF
C
         ELSEIF(ICASE .LT. 0) THEN
            IF    (KT .EQ. 2) THEN
               IESUBR   = 'FFTCAL2I'
               CALL FFTCAL2I(RESULT,STORAG,PHASER,PHASEI,MDIM,M0,NDIM,N0
     &                                           ,ND1,ND2,ND3,NST,MSIGN
     &                                           ,IESUB,IERPLC,IERROR)
            ELSEIF(KT .EQ. 3) THEN
               IESUBR   = 'FFTCAL3I'
               CALL FFTCAL3I(RESULT,STORAG,PHASER,PHASEI,MDIM,M0,NDIM,N0
     &                                           ,ND1,ND2,ND3,NST,MSIGN
     &                                           ,IESUB,IERPLC,IERROR)
            ELSEIF(KT .EQ. 4) THEN
               IESUBR   = 'FFTCAL4I'
               CALL FFTCAL4I(RESULT,STORAG,PHASER,PHASEI,MDIM,M0,NDIM,N0
     &                                           ,ND1,ND2,ND3,NST,MSIGN
     &                                           ,IESUB,IERPLC,IERROR)
            ELSEIF(KT .EQ. 5) THEN
               IESUBR   = 'FFTCAL5I'
               CALL FFTCAL5I(RESULT,STORAG,PHASER,PHASEI,MDIM,M0,NDIM,N0
     &                                           ,ND1,ND2,ND3,NST,MSIGN
     &                                           ,IESUB,IERPLC,IERROR)
            ELSEIF(KT .EQ. 6) THEN
               IESUBR   = 'FFTCAL6I'
               CALL FFTCAL6I(RESULT,STORAG,PHASER,PHASEI,MDIM,M0,NDIM,N0
     &                                           ,ND1,ND2,ND3,NST,MSIGN
     &                                           ,IESUB,IERPLC,IERROR)
            ELSEIF(KT .GT. 6) THEN
               IESUBR   = 'FFTCAL7I'
               CALL FFTCAL7I(RESULT,STORAG,PHASER,PHASEI,MDIM,M0,NDIM,N0
     &                                           ,ND1,ND2,ND3,NST,MSIGN
     &                                           ,IESUB,IERPLC,IERROR)
            ENDIF
         ENDIF
C
         IF(IERROR .NE. 0) THEN
            WRITE(KUNIT,1000) IESUBR,IESUB,IERROR,IERPLC,KT,MSIGN
            RETURN
         ENDIF
C
         MSIGN   = -MSIGN
         IF(ICASE .GT. 0) KPAR    =  KPAR*KT
         IF(ICASE .LT. 0) NPAR    =  NPAR*KT
      ENDDO
C
C
      IF(ICASE .GT. 0) FACTOR = 1.0/SQRT(DFLOAT(NP))
      IF(ICASE .LT. 0) FACTOR = 2.0/SQRT(DFLOAT(NP))
      DO KV      = 1, NP
      DO IV      = 1, MP
      IF(MSIGN .GT. 0) RESULT(IV,KV) = FACTOR*RESULT(IV,KV)
      IF(MSIGN .LT. 0) RESULT(IV,KV) = FACTOR*STORAG(IV,KV)
      ENDDO
      ENDDO
C
C
C
C 4.0 RETURN AND END
C
      RETURN
 1000 FORMAT(/,5X,'FFT ERROR IN ',A8,2x,'IESUB   = ',a8
     &        ,2X,': IERROR   = ',I5
     &        ,2X,'AT IERPLC  = ',I5,7X,'FOR KT      = ',I5
     &        ,2X,'AND MSIGN  = ',I2)
      END
      SUBROUTINE FFTCAL2F(FREAL,FIMAG,PHASER,PHASEI,M0DIM,M0,N0DIM,N0
     &                                      ,ND1,ND2,ND3,NST,SWITCH
     &                                      ,IERSUB,IERPLC,IERROR)
C
C
C ------------------------------------------------------------------------------------
C
C CALCULATE FAST FOURIER TRANSFORM OF FREAL AND STORE IN FIMAG
C COOLEY-TOOKEY ALGORITHM FOR SPECIAL CASE WITH ND2 = 2
C
C ----------------------------------------------------------------------------------------
C
      IMPLICIT NONE
C
      CHARACTER*16 SUBNAM
      INTEGER      I,        IV,       K,        KV,       KPV
      INTEGER      IVP,      KVP,      LVP,      J,        JV
      INTEGER      KVI,      KVO,      LVI,      LVO
      INTEGER      ND3P,     ND3Q,     ND3QP
      INTEGER      NDTL0,    NDTL1,    NDTOTL0,  NDTOTL1
      INTEGER      N3PARITY
      REAL*8       FR1,      FR2
      REAL*8       ZVAR1,    ZVAR2
      REAL*8       ZVAR1R,   ZVAR1I,   ZVAR2R,   ZVAR2I
      REAL*8       CSK1,     CSK2
      REAL*8       SNK1,     SNK2
C
      CHARACTER*16 IERSUB
      INTEGER      M0,       N0
      INTEGER      M0DIM,    N0DIM
      INTEGER      ND1,      ND2,      ND3,      NST,     SWITCH
      INTEGER      IERPLC,    IERROR
C
      REAL*8       FREAL   (M0DIM,N0DIM),FIMAG   (M0DIM,N0DIM)
      REAL*8       PHASER  (N0DIM),      PHASEI  (N0DIM)
C
C
      REAL*8,      DIMENSION(:,:,:),ALLOCATABLE:: WORKINP
      REAL*8,      DIMENSION(:,:,:),ALLOCATABLE:: WORKOUT
      REAL*8,      DIMENSION(:,:),  ALLOCATABLE:: CST
      REAL*8,      DIMENSION(:,:),  ALLOCATABLE:: SNT
C
C
C
C 1.0 INITIALIZATION
C
      SUBNAM    = 'FFTCAL2F'
      IERSUB    = ''
      IERPLC    = 0
      IERROR    = 0
C
      ND3P      = (ND3-1)/2
      ND3Q      =  ND3/2
      ND3QP     =  ND3Q  + 1
      N3PARITY  =  ND3   - 2*(ND3 /2)
C
      IF(ND2 .NE. 2) THEN
         IERSUB   =  SUBNAM
         IERPLC   = -1
         IERROR   =  1
         RETURN
      ENDIF
C
      IF(IABS(SWITCH) .NE. +1) THEN
         IERSUB   =  SUBNAM
         IERPLC   = -2
         IERROR   =  SWITCH
         RETURN
      ENDIF
C
C
C
C 2.0 STORE THE INPUT DATA IN DIFFERENT SIZE WORK ARRAYS
C
C 2.1 CHECK FOR CONSISTENT DIMENSIONS
C
      NDTOTL0   = M0*N0
      NDTOTL1   = ND1*ND2*ND3
      IF(NDTOTL0 .NE. NDTOTL1) THEN
         IERSUB   =  SUBNAM
         IERPLC   = +1
         IERROR   = NDTOTL0 - NDTOTL1
         RETURN
      ENDIF
C
      NDTL0     = ND2*ND3
      NDTL1     = N0 - NST + 1
      IF(NDTL0 .NE. NDTL1) THEN
         IERSUB   =  SUBNAM
         IERPLC   = +2
         IERROR   = NDTL0 - NDTL1
         RETURN
      ENDIF
C
      IF(NST .GT. N0) THEN
         IERSUB   =  SUBNAM
         IERPLC   = +3
         IERROR   =  NST
         RETURN
      ENDIF
C
C
C 2.2 ALLOCATE INTERNAL ARRAYS
C
      ALLOCATE (WORKINP(ND1,ND2,ND3))
      ALLOCATE (WORKOUT(ND1,ND3,ND2))
      ALLOCATE (CST    (ND3,ND2))
      ALLOCATE (SNT    (ND3,ND2))
C
C
C 2.3 FILL IN THE INTERNAL ARRAYS
C
C 2.3.1 FILL IN THE WORK ARRAYS WITH DIMENSIONS (ND1,ND2,ND3) AND (ND1,ND3,ND2)
C
      IVP      = 0
      KVI      = 1
      LVI      = 1
      KVO      = 1
      LVO      = 1
      DO J     = 1,N0
      JV       = J
      DO I     = 1,M0
      IV       = I
      IVP      = IVP + 1
      IF(IVP .GT. ND1) THEN
         IVP      = 1
         KVI      = KVI + 1
         LVO      = LVO + 1
         IF(KVI .GT. ND2) THEN
            KVI       = 1
            LVI       = LVI + 1
            IF(LVI .GT. ND3) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +4
               IERROR   = LVI
               RETURN
            ENDIF
         ENDIF
         IF(LVO .GT. ND3) THEN
            LVO       = 1
            KVO       = KVO + 1
            IF(KVO .GT. ND2) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +5
               IERROR   = KVO
               RETURN
            ENDIF
         ENDIF
      ENDIF
C
      IF    (SWITCH .EQ. +1) THEN
         FR1    = FREAL(IV,JV)
      ELSEIF(SWITCH .EQ. -1) THEN
         FR1    = FIMAG(IV,JV)
      ENDIF
      WORKINP(IVP,KVI,LVI)  = FR1
      WORKOUT(IVP,LVO,KVO)  = 0.0
      ENDDO
      ENDDO
C
C 2.3.2 FILL IN THE PHASE ARRAYS FROM K = NST TO N0
C
      LVP      = 0
      KVP      = 1
      DO J     = NST,N0
      JV       = J
      LVP      = LVP + 1
      IF(LVP .GT. ND3) THEN
         LVP      = 1
         KVP      = KVP + 1
         IF(KVP .GT. ND2) THEN
            IERSUB   =  SUBNAM
            IERPLC   = +6
            IERROR   = KVP
            RETURN
         ENDIF
      ENDIF
      CST(LVP,KVP)  = PHASER(JV)
      SNT(LVP,KVP)  = PHASEI(JV)
      ENDDO
C
C
C
C 3.0 KV = 1: SET FIRST ELEMENT OF MIDDLE INDEX
C
      DO IV              = 1,ND1
      ZVAR1              =           WORKINP(IV,1,1)  +
     &                               WORKINP(IV,2,1)
      ZVAR2              =           WORKINP(IV,1,1)  -
     &                               WORKINP(IV,2,1)
      WORKOUT(IV,1,1)    =          +ZVAR1
      WORKOUT(IV,1,2)    =          +ZVAR2
      ENDDO
C
C
C
C 4.0 FOR KV = 2, ND3P + 1
C
C 4.1 FOR ND1 > (ND3-1)/2
C
      IF    (ND1 .GT. ND3P) THEN
         DO IV              = 1,ND1
            DO K               = 1,ND3P
            KV                 = K   + 1
            KPV                = ND3 - K + 1
            CSK1               =  1.0
            CSK2               =  CST(KV,1)
            SNK1               =  0.0
            SNK2               =  SNT(KV,1)
C
            ZVAR1R             =  CSK1*WORKINP(IV,1, KV)  -
     &                                           SNK1*WORKINP(IV,1, KPV)
            ZVAR1I             =  CSK1*WORKINP(IV,1, KPV) +
     &                                           SNK1*WORKINP(IV,1, KV)
            ZVAR2R             =  CSK2*WORKINP(IV,2, KV)  -
     &                                           SNK2*WORKINP(IV,2, KPV)
            ZVAR2I             =  CSK2*WORKINP(IV,2, KPV) +
     &                                           SNK2*WORKINP(IV,2, KV)
C
            WORKOUT(IV,KV, 1)  =      ZVAR1R   +      ZVAR2R
            WORKOUT(IV,KPV,1)  =      ZVAR1R   -      ZVAR2R
            WORKOUT(IV,KV, 2)  =     -ZVAR1I   +      ZVAR2I
            WORKOUT(IV,KPV,2)  =      ZVAR1I   +      ZVAR2I
            ENDDO
         ENDDO
C
C
C 4.2 FOR ND1 > (ND3-1)/2
C
      ELSEIF(ND1 .LE. ND3P) THEN
         DO K               = 1,ND3P
         KV                 = K   + 1
         KPV                = ND3 - K + 1
         CSK1               = 1.0
         CSK2               = CST(KV,1)
         SNK1               = 0.0
         SNK2               = SNT(KV,1)
            DO IV              = 1,ND1
C
            ZVAR1R             =  CSK1*WORKINP(IV,1, KV)  -
     &                                           SNK1*WORKINP(IV,1, KPV)
            ZVAR1I             =  CSK1*WORKINP(IV,1, KPV) +
     &                                           SNK1*WORKINP(IV,1, KV)
            ZVAR2R             =  CSK2*WORKINP(IV,2, KV)  -
     &                                           SNK2*WORKINP(IV,2, KPV)
            ZVAR2I             =  CSK2*WORKINP(IV,2, KPV) +
     &                                           SNK2*WORKINP(IV,2, KV)
C
            WORKOUT(IV,KV, 1)  =      ZVAR1R   +      ZVAR2R
            WORKOUT(IV,KPV,1)  =      ZVAR1R   -      ZVAR2R
            WORKOUT(IV,KV, 2)  =     -ZVAR1I   +      ZVAR2I
            WORKOUT(IV,KPV,2)  =      ZVAR1I   +      ZVAR2I
            ENDDO
         ENDDO
      END IF
C
C
C 4.3 SPECIAL CASE FOR K = ND3/2 WITH ND3 EVEN SHOULD NOT BE NEEDED
C
      IF(N3PARITY .EQ. 0) THEN
         IERSUB   =  SUBNAM
         IERPLC   = -3
         IERROR   =  4
C
         DO IV                 = 1,ND1
         WORKOUT(IV,ND3QP,1 )  = +WORKINP(IV,1, ND3QP)
         WORKOUT(IV,ND3QP,2 )  = -WORKINP(IV,2, ND3QP)
         ENDDO
      ENDIF
C
C
C
C 5.0 SAVE THE RESULTS AND CLEAN UP
C
C 5.1 SAVE THE RESULTS
C
      IVP      = 0
      KVI      = 1
      LVI      = 1
      KVO      = 1
      LVO      = 1
      DO J     = 1,N0
      JV       = J
      DO I     = 1,M0
      IV       = I
      IVP      = IVP + 1
      IF(IVP .GT. ND1) THEN
         IVP      = 1
         KVI      = KVI + 1
         LVO      = LVO + 1
         IF(KVI .GT. ND2) THEN
            KVI       = 1
            LVI       = LVI + 1
            IF(LVI .GT. ND3) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +7
               IERROR   = LVI
               RETURN
            ENDIF
         ENDIF
         IF(LVO .GT. ND3) THEN
            LVO       = 1
            KVO       = KVO + 1
            IF(KVO .GT. ND2) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +8
               IERROR   = KVO
               RETURN
            ENDIF
         ENDIF
      ENDIF
C
      FR1    = WORKINP(IVP,KVI,LVI)
      FR2    = WORKOUT(IVP,LVO,KVO)
      IF    (SWITCH .EQ. +1) THEN
         FREAL(IV,JV)  = FR1
         FIMAG(IV,JV)  = FR2
      ELSEIF(SWITCH .EQ. -1) THEN
         FREAL(IV,JV)  = FR2
         FIMAG(IV,JV)  = FR1
      ENDIF
      ENDDO
      ENDDO
C
C
C 5.2 DEALLOCATE INTERNAL ARRAYS
C
      DEALLOCATE(WORKINP)
      DEALLOCATE(WORKOUT)
      DEALLOCATE(CST)
      DEALLOCATE(SNT)
C
C
C
C 6.0 RETURN AND END
C
      RETURN
      END
      SUBROUTINE FFTCAL2I(FREAL,FIMAG,PHASER,PHASEI,M0DIM,M0,N0DIM,N0
     &                                      ,ND1,ND2,ND3,NST,SWITCH
     &                                      ,IERSUB,IERPLC,IERROR)
C
C
C ------------------------------------------------------------------------------------
C
C CALCULATE INVERSE FAST FOURIER TRANSFORM OF FREAL AND STORE IN FIMAG
C COOLEY-TOOKEY ALGORITHM FOR SPECIAL CASE WITH ND2 = 2
C
C ----------------------------------------------------------------------------------------
C
      IMPLICIT NONE
C
      CHARACTER*16 SUBNAM
      INTEGER      I,        IV,       K,        KV,       KPV
      INTEGER      IVP,      KVP,      LVP,      J,        JV
      INTEGER      KVI,      KVO,      LVI,      LVO
      INTEGER      ND3P,     ND3Q,     ND3QP
      INTEGER      NDTL0,    NDTL1,    NDTOTL0,  NDTOTL1
      INTEGER      N3PARITY
      REAL*8       FR1,      FR2
      REAL*8       ZVAR1,    ZVAR2
      REAL*8       ZVAR1R,   ZVAR1I,   ZVAR2R,   ZVAR2I
      REAL*8       CSK1,     CSK2
      REAL*8       SNK1,     SNK2
C
      CHARACTER*16 IERSUB
      INTEGER      M0,       N0
      INTEGER      M0DIM,    N0DIM
      INTEGER      ND1,      ND2,      ND3,      NST,      SWITCH
      INTEGER      IERPLC,   IERROR
C
      REAL*8       FREAL   (M0DIM,N0DIM),FIMAG   (M0DIM,N0DIM)
      REAL*8       PHASER  (N0DIM),      PHASEI  (N0DIM)
C
C
      REAL*8,      DIMENSION(:,:,:),ALLOCATABLE:: WORKINP
      REAL*8,      DIMENSION(:,:,:),ALLOCATABLE:: WORKOUT
      REAL*8,      DIMENSION(:,:),  ALLOCATABLE:: CST
      REAL*8,      DIMENSION(:,:),  ALLOCATABLE:: SNT
C
C
C 1.0 INITIALIZATION
C
      SUBNAM    = 'FFTCAL2I'
      IERSUB    = ''
      IERPLC    = 0
      IERROR    = 0
C
      ND3P      = (ND3-1)/2
      ND3Q      =  ND3/2
      ND3QP     =  ND3Q  + 1
      N3PARITY  =  ND3   - 2*(ND3 /2)
C
      IF(ND2 .NE. 2) THEN
         IERSUB   =  SUBNAM
         IERPLC   = -1
         IERROR   =  1
         RETURN
      ENDIF
C
      IF(IABS(SWITCH) .NE. +1) THEN
         IERSUB   =  SUBNAM
         IERPLC   = -2
         IERROR   =  SWITCH
         RETURN
      ENDIF
C
C
C
C 2.0 STORE THE INPUT DATA IN DIFFERENT SIZE WORK ARRAYS
C
C 2.1 CHECK FOR CONSISTENT DIMENSIONS
C
      NDTOTL0   = M0*N0
      NDTOTL1   = ND1*ND2*ND3
      IF(NDTOTL0 .NE. NDTOTL1) THEN
         IERSUB   =  SUBNAM
         IERPLC   = +1
         IERROR   = NDTOTL0 - NDTOTL1
         RETURN
      ENDIF
C
      NDTL0     = ND2*ND3
      NDTL1     = N0 - NST + 1
      IF(NDTL0 .NE. NDTL1) THEN
         IERSUB   =  SUBNAM
         IERPLC   = +2
         IERROR   = NDTL0 - NDTL1
         RETURN
      ENDIF
C
      IF(NST .GT. N0) THEN
         IERSUB   =  SUBNAM
         IERPLC   = +3
         IERROR   =  NST
         RETURN
      ENDIF
C
C
C 2.2 ALLOCATE INTERNAL ARRAYS
C
      ALLOCATE (WORKINP(ND1,ND3,ND2))
      ALLOCATE (WORKOUT(ND1,ND2,ND3))
      ALLOCATE (CST    (ND3,ND2))
      ALLOCATE (SNT    (ND3,ND2))
C
C
C 2.3 FILL IN THE INTERNAL ARRAYS
C
C 2.3.1 FILL IN THE WORK ARRAYS WITH DIMENSIONS (ND1,ND3,ND2) AND (ND1,ND2,ND3)
C
      IVP      = 0
      KVI      = 1
      LVI      = 1
      KVO      = 1
      LVO      = 1
      DO J     = 1,N0
      JV       = J
      DO I     = 1,M0
      IV       = I
      IVP      = IVP + 1
      IF(IVP .GT. ND1) THEN
         IVP      = 1
         LVI      = LVI + 1
         KVO      = KVO + 1
         IF(LVI .GT. ND3) THEN
            LVI       = 1
            KVI       = KVI + 1
            IF(KVI .GT. ND2) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +4
               IERROR   = KVI
               RETURN
            ENDIF
         ENDIF
         IF(KVO .GT. ND2) THEN
            KVO       = 1
            LVO       = LVO + 1
            IF(LVO .GT. ND3) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +5
               IERROR   = LVO
               RETURN
            ENDIF
         ENDIF
      ENDIF
C
      IF    (SWITCH .EQ. +1) THEN
         FR1    = FREAL(IV,JV)
      ELSEIF(SWITCH .EQ. -1) THEN
         FR1    = FIMAG(IV,JV)
      ENDIF
      WORKINP(IVP,LVI,KVI)  = FR1
      WORKOUT(IVP,KVO,LVO)  = 0.0
      ENDDO
      ENDDO
C
C 2.3.2 FILL IN THE PHASE ARRAYS FROM K = NST TO N0
C
      LVP      = 0
      KVP      = 1
      DO J     = NST,N0
      JV       = J
      LVP      = LVP + 1
      IF(LVP .GT. ND3) THEN
         LVP      = 1
         KVP      = KVP + 1
         IF(KVP .GT. ND2) THEN
            IERSUB   =  SUBNAM
            IERPLC   = +6
            IERROR   = KVP
            RETURN
         ENDIF
      ENDIF
      CST(LVP,KVP)  = PHASER(JV)
      SNT(LVP,KVP)  = PHASEI(JV)
      ENDDO
C
C
C
C 3.0 KV = 1: SET FIRST ELEMENT OF MIDDLE INDEX
C
      DO IV              = 1,ND1
      ZVAR1              =           WORKINP(IV,1,1)  +
     &                               WORKINP(IV,1,2)
      ZVAR2              =           WORKINP(IV,1,1)  -
     &                               WORKINP(IV,1,2)
      WORKOUT(IV,1,1)    =          +ZVAR1
      WORKOUT(IV,2,1)    =          +ZVAR2
      ENDDO
C
C
C
C 4.0 FOR KV = 2, ND3P + 1
C
C 4.1 FOR ND1 > (ND3-1)/2
C
      IF    (ND1 .GT. ND3P) THEN
         DO IV              = 1,ND1
            DO K               = 1,ND3P
            KV                 = K   + 1
            KPV                = ND3 - K + 1
            CSK1               = 1.0
            CSK2               = CST(KV,1)
            SNK1               = 0.0
            SNK2               = SNT(KV,1)
C
            ZVAR1R             = +WORKINP(IV,KV,  1)   +
     &                            WORKINP(IV,KPV, 1)
            ZVAR1I             = -WORKINP(IV,KV,  2)   +
     &                            WORKINP(IV,KPV, 2)
            ZVAR2R             = +WORKINP(IV,KV,  1)   -
     &                            WORKINP(IV,KPV, 1)
            ZVAR2I             = +WORKINP(IV,KV,  2)   +
     &                            WORKINP(IV,KPV, 2)
C
            WORKOUT(IV,1, KV)  =    +CSK1*ZVAR1R   -     SNK1*ZVAR1I
            WORKOUT(IV,1, KPV) =    +CSK1*ZVAR1I   +     SNK1*ZVAR1R
            WORKOUT(IV,2, KV)  =    +CSK2*ZVAR2R   -     SNK2*ZVAR2I
            WORKOUT(IV,2, KPV) =    +CSK2*ZVAR2I   +     SNK2*ZVAR2R
            ENDDO
         ENDDO
C
C
C 4.2 FOR ND1 > (ND3-1)/2
C
      ELSEIF(ND1 .LE. ND3P) THEN
         DO K               = 1,ND3P
         KV                 = K   + 1
         KPV                = ND3 - K + 1
         CSK1               = 1.0
         CSK2               = CST(KV,1)
         SNK1               = 0.0
         SNK2               = SNT(KV,1)
            DO IV              = 1,ND1
            ZVAR1R             = +WORKINP(IV,KV,  1)   +
     &                            WORKINP(IV,KPV, 1)
            ZVAR1I             = -WORKINP(IV,KV,  2)   +
     &                            WORKINP(IV,KPV, 2)
            ZVAR2R             = +WORKINP(IV,KV,  1)   -
     &                            WORKINP(IV,KPV, 1)
            ZVAR2I             = +WORKINP(IV,KV,  2)   +
     &                            WORKINP(IV,KPV, 2)
C
            WORKOUT(IV,1, KV)  =    +CSK1*ZVAR1R   -     SNK1*ZVAR1I
            WORKOUT(IV,1, KPV) =    +CSK1*ZVAR1I   +     SNK1*ZVAR1R
            WORKOUT(IV,2, KV)  =    +CSK2*ZVAR2R   -     SNK2*ZVAR2I
            WORKOUT(IV,2, KPV) =    +CSK2*ZVAR2I   +     SNK2*ZVAR2R
            ENDDO
         ENDDO
      END IF
C
C
C 4.3 SPECIAL CASE FOR K = ND3/2 WITH ND3 EVEN SHOULD NOT BE NEEDED
C
      IF(N3PARITY .EQ. 0) THEN
         IERSUB   =  SUBNAM
         IERPLC   = -3
         IERROR   =  4
C
         DO IV                 = 1,ND1
         WORKOUT(IV,1, ND3QP)  = +WORKINP(IV,ND3QP,1 )
         WORKOUT(IV,2, ND3QP)  = +WORKINP(IV,ND3QP,2 )
         ENDDO
      ENDIF
C
C
C
C 5.0 SAVE THE RESULTS AND CLEAN UP
C
C 5.1 SAVE THE RESULTS
C
      IVP      = 0
      KVI      = 1
      LVI      = 1
      KVO      = 1
      LVO      = 1
      DO J     = 1,N0
      JV       = J
      DO I     = 1,M0
      IV       = I
      IVP      = IVP + 1
      IF(IVP .GT. ND1) THEN
         IVP      = 1
         LVI      = LVI + 1
         KVO      = KVO + 1
         IF(LVI .GT. ND3) THEN
            LVI       = 1
            KVI       = KVI + 1
            IF(KVI .GT. ND2) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +7
               IERROR   = KVI
               RETURN
            ENDIF
         ENDIF
         IF(KVO .GT. ND2) THEN
            KVO       = 1
            LVO       = LVO + 1
            IF(LVO .GT. ND3) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +8
               IERROR   = LVO
               RETURN
            ENDIF
         ENDIF
      ENDIF
C
      FR1    = WORKINP(IVP,LVI,KVI)
      FR2    = WORKOUT(IVP,KVO,LVO)
      IF    (SWITCH .EQ. +1) THEN
         FREAL(IV,JV)  = FR1
         FIMAG(IV,JV)  = FR2
      ELSEIF(SWITCH .EQ. -1) THEN
         FREAL(IV,JV)  = FR2
         FIMAG(IV,JV)  = FR1
      ENDIF
      ENDDO
      ENDDO
C
C
C 5.2 DEALLOCATE INTERNAL ARRAYS
C
      DEALLOCATE(WORKINP)
      DEALLOCATE(WORKOUT)
      DEALLOCATE(CST)
      DEALLOCATE(SNT)
C
C
C
C 6.0 RETURN AND END
C
      RETURN
      END
      SUBROUTINE FFTCAL3F(FREAL,FIMAG,PHASER,PHASEI,M0DIM,M0,N0DIM,N0
     &                                      ,ND1,ND2,ND3,NST,SWITCH
     &                                      ,IERSUB,IERPLC,IERROR)
C
C
C ------------------------------------------------------------------------------------
C
C CALCULATE FAST FOURIER TRANSFORM OF FREAL AND STORE IN FIMAG
C COOLEY-TOOKEY ALGORITHM FOR SPECIAL CASE WITH ND2 = 3
C
C ----------------------------------------------------------------------------------------
C
      IMPLICIT NONE
C
      CHARACTER*16 SUBNAM
      INTEGER      I,        IV,       K,        KV,       KPV
      INTEGER      IVP,      KVP,      LVP,      J,        JV
      INTEGER      KVI,      KVO,      LVI,      LVO
      INTEGER      ND3P,     ND3Q,     ND3QP
      INTEGER      NDTL0,    NDTL1,    NDTOTL0,  NDTOTL1
      INTEGER      N3PARITY
      REAL*8       ROOT3D2
      REAL*8       FR1,      FR2
      REAL*8       ZVAR0,    ZVAR1,    ZVAR2,    ZVAR3
      REAL*8       ZVAR1R,   ZVAR1I,   ZVAR2R,   ZVAR2I,
     &             ZVAR3R,   ZVAR3I
      REAL*8       ZVAR0R1,  ZVAR0I1,  ZVAR0R2,  ZVAR0I2,
     &             ZVAR0R3,  ZVAR0I3
      REAL*8       CSK1,     CSK2,     CSK3
      REAL*8       SNK1,     SNK2,     SNK3
C
      CHARACTER*16 IERSUB
      INTEGER      M0,       N0
      INTEGER      M0DIM,    N0DIM
      INTEGER      ND1,      ND2,      ND3,      NST,      SWITCH
      INTEGER      IERPLC,   IERROR
C
      REAL*8       FREAL   (M0DIM,N0DIM),FIMAG   (M0DIM,N0DIM)
      REAL*8       PHASER  (N0DIM),      PHASEI  (N0DIM)
C
C
      REAL*8,      DIMENSION(:,:,:),ALLOCATABLE:: WORKINP
      REAL*8,      DIMENSION(:,:,:),ALLOCATABLE:: WORKOUT
      REAL*8,      DIMENSION(:,:),  ALLOCATABLE:: CST
      REAL*8,      DIMENSION(:,:),  ALLOCATABLE:: SNT
C
C
C
C 1.0 INITIALIZATION
C
      SUBNAM    = 'FFTCAL3F'
      IERPLC    = 0
      IERROR    = 0
C
      ND3P      = (ND3-1)/2
      ND3Q      =  ND3/2
      ND3QP     =  ND3Q  + 1
      N3PARITY  =  ND3   - 2*(ND3 /2)
C
      IF(ND2 .NE. 3) THEN
         IERSUB   =  SUBNAM
         IERPLC   = -1
         IERROR   =  1
         RETURN
      ENDIF
C
      IF(IABS(SWITCH) .NE. +1) THEN
         IERSUB   =  SUBNAM
         IERPLC   = -2
         IERROR   =  SWITCH
         RETURN
      ENDIF
C
      ROOT3D2   = 0.866025403784438646763723170752936
C
C
C
C 2.0 STORE THE INPUT DATA IN DIFFERENT SIZE WORK ARRAYS
C
C 2.1 CHECK FOR CONSISTENT DIMENSIONS
C
      NDTOTL0   = M0*N0
      NDTOTL1   = ND1*ND2*ND3
      IF(NDTOTL0 .NE. NDTOTL1) THEN
         IERSUB   =  SUBNAM
         IERPLC   = +1
         IERROR   = NDTOTL0 - NDTOTL1
         RETURN
      ENDIF
C
      NDTL0     = ND2*ND3
      NDTL1     = N0 - NST + 1
      IF(NDTL0 .NE. NDTL1) THEN
         IERSUB   =  SUBNAM
         IERPLC   = +2
         IERROR   = NDTL0 - NDTL1
         RETURN
      ENDIF
C
      IF(NST .GT. N0) THEN
         IERSUB   =  SUBNAM
         IERPLC   = +3
         IERROR   =  NST
         RETURN
      ENDIF
C
C
C 2.2 ALLOCATE INTERNAL ARRAYS
C
      ALLOCATE (WORKINP(ND1,ND2,ND3))
      ALLOCATE (WORKOUT(ND1,ND3,ND2))
      ALLOCATE (CST    (ND3,ND2))
      ALLOCATE (SNT    (ND3,ND2))
C
C
C 2.3 FILL IN THE INTERNAL ARRAYS
C
C 2.3.1 FILL IN THE WORK ARRAYS WITH DIMENSIONS (ND1,ND2,ND3) AND (ND1,ND3,ND2)
C
      IVP      = 0
      KVI      = 1
      LVI      = 1
      KVO      = 1
      LVO      = 1
      DO J     = 1,N0
      JV       = J
      DO I     = 1,M0
      IV       = I
      IVP      = IVP + 1
      IF(IVP .GT. ND1) THEN
         IVP      = 1
         KVI      = KVI + 1
         LVO      = LVO + 1
         IF(KVI .GT. ND2) THEN
            KVI       = 1
            LVI       = LVI + 1
            IF(LVI .GT. ND3) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +4
               IERROR   = LVI
               RETURN
            ENDIF
         ENDIF
         IF(LVO .GT. ND3) THEN
            LVO       = 1
            KVO       = KVO + 1
            IF(KVO .GT. ND2) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +5
               IERROR   = KVO
               RETURN
            ENDIF
         ENDIF
      ENDIF
C
      IF    (SWITCH .EQ. +1) THEN
         FR1    = FREAL(IV,JV)
      ELSEIF(SWITCH .EQ. -1) THEN
         FR1    = FIMAG(IV,JV)
      ENDIF
      WORKINP(IVP,KVI,LVI)  = FR1
      WORKOUT(IVP,LVO,KVO)  = 0.0
      ENDDO
      ENDDO
C
C 2.3.2 FILL IN THE PHASE ARRAYS FROM K = NST TO N0
C
      LVP      = 0
      KVP      = 1
      DO J     = NST,N0
      JV       = J
      LVP      = LVP + 1
      IF(LVP .GT. ND3) THEN
         LVP      = 1
         KVP      = KVP + 1
         IF(KVP .GT. ND2) THEN
            IERSUB   =  SUBNAM
            IERPLC   = +6
            IERROR   = KVP
            RETURN
         ENDIF
      ENDIF
      CST(LVP,KVP)  = PHASER(JV)
      SNT(LVP,KVP)  = PHASEI(JV)
      ENDDO
C
C
C
C 3.0 KV = 1: SET FIRST ELEMENT OF MIDDLE INDEX
C
      DO IV              = 1,ND1
      ZVAR0              =           WORKINP(IV,1,1)
      ZVAR1              =           WORKINP(IV,2,1)  +
     &                               WORKINP(IV,3,1)
      ZVAR2              =           WORKINP(IV,2,1)  -
     &                               WORKINP(IV,3,1)
C
      WORKOUT(IV,1,1)    =          +ZVAR0            +       ZVAR1
      WORKOUT(IV,1,2)    =          +ZVAR0            -   0.5*ZVAR1
      WORKOUT(IV,1,3)    =  -ROOT3D2*ZVAR2
      ENDDO
C
C
C
C 4.0 FOR KV = 2, ND3P + 1
C
C 4.1 FOR ND1 > (ND3-1)/2
C
      IF    (ND1 .GT. ND3P) THEN
         DO IV              = 1,ND1
            DO K               = 1,ND3P
            KV                 = K   + 1
            KPV                = ND3 - K + 1
            CSK1               =  1.0
            CSK2               =  CST(KV,1)
            CSK3               =  CST(KV,2)
            SNK1               =  0.0
            SNK2               =  SNT(KV,1)
            SNK3               =  SNT(KV,2)
C
            ZVAR1R             =  CSK1*WORKINP(IV,1, KV)  -
     &                                           SNK1*WORKINP(IV,1, KPV)
            ZVAR1I             =  CSK1*WORKINP(IV,1, KPV) +
     &                                           SNK1*WORKINP(IV,1, KV)
            ZVAR2R             =  CSK2*WORKINP(IV,2, KV)  -
     &                                           SNK2*WORKINP(IV,2, KPV)
            ZVAR2I             =  CSK2*WORKINP(IV,2, KPV) +
     &                                           SNK2*WORKINP(IV,2, KV)
            ZVAR3R             =  CSK3*WORKINP(IV,3, KV)  -
     &                                           SNK3*WORKINP(IV,3, KPV)
            ZVAR3I             =  CSK3*WORKINP(IV,3, KPV) +
     &                                           SNK3*WORKINP(IV,3, KV)
C
            ZVAR0R1            =           ZVAR2R   +      ZVAR3R
            ZVAR0I1            =           ZVAR2I   +      ZVAR3I
            ZVAR0R2            =           ZVAR1R   -  0.5*ZVAR0R1
            ZVAR0I2            =           ZVAR1I   -  0.5*ZVAR0I1
            ZVAR0R3            = +ROOT3D2*(ZVAR2R   -      ZVAR3R)
            ZVAR0I3            = +ROOT3D2*(ZVAR2I   -      ZVAR3I)
C
            WORKOUT(IV,KV, 1)  =      ZVAR1R   +      ZVAR0R1
            WORKOUT(IV,KPV,1)  =      ZVAR0R2  -      ZVAR0I3
            WORKOUT(IV,KV, 2)  =      ZVAR0R2  +      ZVAR0I3
            WORKOUT(IV,KPV,2)  =      ZVAR0I2  -      ZVAR0R3
            WORKOUT(IV,KV, 3)  =     -ZVAR0I2  -      ZVAR0R3
            WORKOUT(IV,KPV,3)  =      ZVAR1I   +      ZVAR0I1
            ENDDO
         ENDDO
C
C
C 4.2 FOR ND1 > (ND3-1)/2
C
      ELSEIF(ND1 .LE. ND3P) THEN
         DO K               = 1,ND3P
         KV                 = K   + 1
         KPV                = ND3 - K + 1
         CSK1               = 1.0
         CSK2               = CST(KV,1)
         CSK3               = CST(KV,2)
         SNK1               = 0.0
         SNK2               = SNT(KV,1)
         SNK3               = SNT(KV,2)
            DO IV              = 1,ND1
            ZVAR1R             =  CSK1*WORKINP(IV,1, KV)  -
     &                                           SNK1*WORKINP(IV,1, KPV)
            ZVAR1I             =  CSK1*WORKINP(IV,1, KPV) +
     &                                           SNK1*WORKINP(IV,1, KV)
            ZVAR2R             =  CSK2*WORKINP(IV,2, KV)  -
     &                                           SNK2*WORKINP(IV,2, KPV)
            ZVAR2I             =  CSK2*WORKINP(IV,2, KPV) +
     &                                           SNK2*WORKINP(IV,2, KV)
            ZVAR3R             =  CSK3*WORKINP(IV,3, KV)  -
     &                                           SNK3*WORKINP(IV,3, KPV)
            ZVAR3I             =  CSK3*WORKINP(IV,3, KPV) +
     &                                           SNK3*WORKINP(IV,3, KV)
C
            ZVAR0R1            =           ZVAR2R   +      ZVAR3R
            ZVAR0I1            =           ZVAR2I   +      ZVAR3I
            ZVAR0R2            =           ZVAR1R   -  0.5*ZVAR0R1
            ZVAR0I2            =           ZVAR1I   -  0.5*ZVAR0I1
            ZVAR0R3            = +ROOT3D2*(ZVAR2R   -      ZVAR3R)
            ZVAR0I3            = +ROOT3D2*(ZVAR2I   -      ZVAR3I)
C
            WORKOUT(IV,KV, 1)  =      ZVAR1R   +      ZVAR0R1
            WORKOUT(IV,KPV,1)  =      ZVAR0R2  -      ZVAR0I3
            WORKOUT(IV,KV, 2)  =      ZVAR0R2  +      ZVAR0I3
            WORKOUT(IV,KPV,2)  =      ZVAR0I2  -      ZVAR0R3
            WORKOUT(IV,KV, 3)  =     -ZVAR0I2  -      ZVAR0R3
            WORKOUT(IV,KPV,3)  =      ZVAR1I   +      ZVAR0I1
            ENDDO
         ENDDO
      END IF
C
C
C 4.3 SPECIAL CASE FOR K = ND3/2 WITH ND3 EVEN SHOULD NOT BE NEEDED
C
      IF(N3PARITY .EQ. 0) THEN
         IERSUB   =  SUBNAM
         IERPLC   = -3
         IERROR   =  4
C
         DO IV                 = 1,ND1
         ZVAR0                 =           WORKINP(IV,1,ND3QP)
         ZVAR1                 =           WORKINP(IV,2,ND3QP) -
     &                                     WORKINP(IV,3,ND3QP)
         ZVAR2                 =           WORKINP(IV,2,ND3QP)
         ZVAR3                 = -ROOT3D2*(WORKINP(IV,2,ND3QP) +
     &                                     WORKINP(IV,3,ND3QP))
C
         WORKOUT(IV,ND3QP, 1)  =          +ZVAR0    + 0.5*ZVAR1
         WORKOUT(IV,ND3QP, 2)  =          +ZVAR0    -     ZVAR1
         WORKOUT(IV,ND3QP, 3)  =          +ZVAR3
         ENDDO
      ENDIF
C
C
C
C 5.0 SAVE THE RESULTS AND CLEAN UP
C
C 5.1 SAVE THE RESULTS
C
      IVP      = 0
      KVI      = 1
      LVI      = 1
      KVO      = 1
      LVO      = 1
      DO J     = 1,N0
      JV       = J
      DO I     = 1,M0
      IV       = I
      IVP      = IVP + 1
      IF(IVP .GT. ND1) THEN
         IVP      = 1
         KVI      = KVI + 1
         LVO      = LVO + 1
         IF(KVI .GT. ND2) THEN
            KVI       = 1
            LVI       = LVI + 1
            IF(LVI .GT. ND3) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +7
               IERROR   = LVI
               RETURN
            ENDIF
         ENDIF
         IF(LVO .GT. ND3) THEN
            LVO       = 1
            KVO       = KVO + 1
            IF(KVO .GT. ND2) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +8
               IERROR   = KVO
               RETURN
            ENDIF
         ENDIF
      ENDIF
C
      FR1    = WORKINP(IVP,KVI,LVI)
      FR2    = WORKOUT(IVP,LVO,KVO)
      IF    (SWITCH .EQ. +1) THEN
         FREAL(IV,JV)  = FR1
         FIMAG(IV,JV)  = FR2
      ELSEIF(SWITCH .EQ. -1) THEN
         FREAL(IV,JV)  = FR2
         FIMAG(IV,JV)  = FR1
      ENDIF
      ENDDO
      ENDDO
C
C
C 5.2 DEALLOCATE INTERNAL ARRAYS
C
      DEALLOCATE(WORKINP)
      DEALLOCATE(WORKOUT)
      DEALLOCATE(CST)
      DEALLOCATE(SNT)
C
C
C
C 6.0 RETURN AND END
C
      RETURN
      END
      SUBROUTINE FFTCAL3I(FREAL,FIMAG,PHASER,PHASEI,M0DIM,M0,N0DIM,N0
     &                                      ,ND1,ND2,ND3,NST,SWITCH
     &                                      ,IERSUB,IERPLC,IERROR)
C
C
C ------------------------------------------------------------------------------------
C
C CALCULATE INVERSE FAST FOURIER TRANSFORM OF FREAL AND STORE IN FIMAG
C COOLEY-TOOKEY ALGORITHM FOR SPECIAL CASE WITH ND2 = 3
C
C ----------------------------------------------------------------------------------------
C
      IMPLICIT NONE
C
      CHARACTER*16 SUBNAM
      INTEGER      I,        IV,       K,        KV,       KPV
      INTEGER      IVP,      KVP,      LVP,      J,        JV
      INTEGER      KVI,      KVO,      LVI,      LVO
      INTEGER      ND3P,     ND3Q,     ND3QP
      INTEGER      NDTL0,    NDTL1,    NDTOTL0,  NDTOTL1
      INTEGER      N3PARITY
      REAL*8       ROOT3D2
      REAL*8       FR1,      FR2
      REAL*8       ZVAR0,    ZVAR1,    ZVAR2,    ZVAR3,    ZVAR4
      REAL*8       ZVAR1R,   ZVAR1I,   ZVAR2R,   ZVAR2I,
     &             ZVAR3R,   ZVAR3I
      REAL*8       ZVAR0R1,  ZVAR0I1,  ZVAR0R2,  ZVAR0I2,
     &             ZVAR0R3,  ZVAR0I3
      REAL*8       ZVAR1R1,  ZVAR1I1,  ZVAR1R2,  ZVAR1I2,
     &             ZVAR1R3,  ZVAR1I3
      REAL*8       CSK1,     CSK2,     CSK3
      REAL*8       SNK1,     SNK2,     SNK3
C
      CHARACTER*16 IERSUB
      INTEGER      M0,       N0
      INTEGER      M0DIM,    N0DIM
      INTEGER      ND1,      ND2,      ND3,      NST,     SWITCH
      INTEGER      IERPLC,   IERROR
C
      REAL*8       FREAL   (M0DIM,N0DIM),FIMAG   (M0DIM,N0DIM)
      REAL*8       PHASER  (N0DIM),      PHASEI  (N0DIM)
C
C
      REAL*8,      DIMENSION(:,:,:),ALLOCATABLE:: WORKINP
      REAL*8,      DIMENSION(:,:,:),ALLOCATABLE:: WORKOUT
      REAL*8,      DIMENSION(:,:),  ALLOCATABLE:: CST
      REAL*8,      DIMENSION(:,:),  ALLOCATABLE:: SNT
C
C
C
C 1.0 INITIALIZATION
C
      SUBNAM    = 'FFTCAL3F'
      IERSUB    = ''
      IERPLC    = 0
      IERROR    = 0
C
      ND3P      = (ND3-1)/2
      ND3Q      =  ND3/2
      ND3QP     =  ND3Q  + 1
      N3PARITY  =  ND3   - 2*(ND3 /2)
C
      IF(ND2 .NE. 3) THEN
         IERSUB   =  SUBNAM
         IERPLC   = -1
         IERROR   =  1
         RETURN
      ENDIF
C
      IF(IABS(SWITCH) .NE. +1) THEN
         IERSUB   =  SUBNAM
         IERPLC   = -2
         IERROR   =  SWITCH
         RETURN
      ENDIF
C
      ROOT3D2   = 0.866025403784438646763723170752936
C
C
C
C 2.0 STORE THE INPUT DATA IN DIFFERENT SIZE WORK ARRAYS
C
C 2.1 CHECK FOR CONSISTENT DIMENSIONS
C
      NDTOTL0   = M0*N0
      NDTOTL1   = ND1*ND2*ND3
      IF(NDTOTL0 .NE. NDTOTL1) THEN
         IERSUB   =  SUBNAM
         IERPLC   = +1
         IERROR   = NDTOTL0 - NDTOTL1
         RETURN
      ENDIF
C
      NDTL0     = ND2*ND3
      NDTL1     = N0 - NST + 1
      IF(NDTL0 .NE. NDTL1) THEN
         IERSUB   =  SUBNAM
         IERPLC   = +2
         IERROR   = NDTL0 - NDTL1
         RETURN
      ENDIF
C
      IF(NST .GT. N0) THEN
         IERSUB   =  SUBNAM
         IERPLC   = +3
         IERROR   =  NST
         RETURN
      ENDIF
C
C
C 2.2 ALLOCATE INTERNAL ARRAYS
C
      ALLOCATE (WORKINP(ND1,ND3,ND2))
      ALLOCATE (WORKOUT(ND1,ND2,ND3))
      ALLOCATE (CST    (ND3,ND2))
      ALLOCATE (SNT    (ND3,ND2))
C
C
C 2.3 FILL IN THE INTERNAL ARRAYS
C
C 2.3.1 FILL IN THE WORK ARRAYS WITH DIMENSIONS (ND1,ND3,ND2) AND (ND1,ND2,ND3)
C
      IVP      = 0
      KVI      = 1
      LVI      = 1
      KVO      = 1
      LVO      = 1
      DO J     = 1,N0
      JV       = J
      DO I     = 1,M0
      IV       = I
      IVP      = IVP + 1
      IF(IVP .GT. ND1) THEN
         IVP      = 1
         LVI      = LVI + 1
         KVO      = KVO + 1
         IF(LVI .GT. ND3) THEN
            LVI       = 1
            KVI       = KVI + 1
            IF(KVI .GT. ND2) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +4
               IERROR   = KVI
               RETURN
            ENDIF
         ENDIF
         IF(KVO .GT. ND2) THEN
            KVO       = 1
            LVO       = LVO + 1
            IF(LVO .GT. ND3) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +5
               IERROR   = LVO
               RETURN
            ENDIF
         ENDIF
      ENDIF
C
      IF    (SWITCH .EQ. +1) THEN
         FR1    = FREAL(IV,JV)
      ELSEIF(SWITCH .EQ. -1) THEN
         FR1    = FIMAG(IV,JV)
      ENDIF
      WORKINP(IVP,LVI,KVI)  = FR1
      WORKOUT(IVP,KVO,LVO)  = 0.0
      ENDDO
      ENDDO
C
C 2.3.2 FILL IN THE PHASE ARRAYS FROM K = NST TO N0
C
      LVP      = 0
      KVP      = 1
      DO J     = NST,N0
      JV       = J
      LVP      = LVP + 1
      IF(LVP .GT. ND3) THEN
         LVP      = 1
         KVP      = KVP + 1
         IF(KVP .GT. ND2) THEN
            IERSUB   =  SUBNAM
            IERPLC   = +6
            IERROR   = KVP
            RETURN
         ENDIF
      ENDIF
      CST(LVP,KVP)  = PHASER(JV)
      SNT(LVP,KVP)  = PHASEI(JV)
      ENDDO
C
C
C
C 3.0 KV = 1: SET FIRST ELEMENT OF MIDDLE INDEX
C
      DO IV              = 1,ND1
      ZVAR0              =           WORKINP(IV,1,1)
      ZVAR1              =           WORKINP(IV,1,2)
      ZVAR2              =           WORKINP(IV,1,3)
      ZVAR3              =           ZVAR0            -   0.5*ZVAR1
      ZVAR4              = +ROOT3D2* WORKINP(IV,1,3)
C
      WORKOUT(IV,1,1)    =          +ZVAR0            +       ZVAR1
      WORKOUT(IV,2,1)    =          +ZVAR3            +       ZVAR4
      WORKOUT(IV,3,1)    =          +ZVAR3            -       ZVAR4
      ENDDO
C
C
C
C 4.0 FOR KV = 2, ND3P + 1
C
C 4.1 FOR ND1 > (ND3-1)/2
C
      IF    (ND1 .GT. ND3P) THEN
         DO IV              = 1,ND1
            DO K               = 1,ND3P
            KV                 = K   + 1
            KPV                = ND3 - K + 1
            CSK1               = 1.0
            CSK2               = CST(KV,1)
            CSK3               = CST(KV,2)
            SNK1               = 0.0
            SNK2               = SNT(KV,1)
            SNK3               = SNT(KV,2)
C
            ZVAR1R             =  WORKINP(IV,KV,  1)
            ZVAR1I             =  WORKINP(IV,KPV, 1)
            ZVAR2R             =  WORKINP(IV,KV,  2)
            ZVAR2I             =  WORKINP(IV,KPV, 2)
            ZVAR3R             =  WORKINP(IV,KV,  3)
            ZVAR3I             =  WORKINP(IV,KPV, 3)
C
            ZVAR0R1            =          +ZVAR2R      +      ZVAR1I
            ZVAR0I1            =          -ZVAR3R      +      ZVAR2I
            ZVAR0R2            =          +ZVAR1R      -  0.5*ZVAR0R1
            ZVAR0I2            =          +ZVAR3I      -  0.5*ZVAR0I1
            ZVAR0R3            = +ROOT3D2*(ZVAR2R      -      ZVAR1I)
            ZVAR0I3            = +ROOT3D2*(ZVAR3R      +      ZVAR2I)
C
            ZVAR1R1            =          +ZVAR1R      +      ZVAR0R1
            ZVAR1I1            =          +ZVAR3I      +      ZVAR0I1
            ZVAR1R2            =          +ZVAR0R2     +      ZVAR0I3
            ZVAR1I2            =          +ZVAR0I2     -      ZVAR0R3
            ZVAR1R3            =          +ZVAR0R2     -      ZVAR0I3
            ZVAR1I3            =          +ZVAR0I2     +      ZVAR0R3
C
            WORKOUT(IV,1, KV)  =    +CSK1*ZVAR1R1   -     SNK1*ZVAR1I1
            WORKOUT(IV,1, KPV) =    +CSK1*ZVAR1I1   +     SNK1*ZVAR1R1
            WORKOUT(IV,2, KV)  =    +CSK2*ZVAR1R2   -     SNK2*ZVAR1I2
            WORKOUT(IV,2, KPV) =    +CSK2*ZVAR1I2   +     SNK2*ZVAR1R2
            WORKOUT(IV,3, KV)  =    +CSK3*ZVAR1R3   -     SNK3*ZVAR1I3
            WORKOUT(IV,3, KPV) =    +CSK3*ZVAR1I3   +     SNK3*ZVAR1R3
            ENDDO
         ENDDO
C
C
C 4.2 FOR ND1 > (ND3-1)/2
C
      ELSEIF(ND1 .LE. ND3P) THEN
         DO K               = 1,ND3P
         KV                 = K   + 1
         KPV                = ND3 - K + 1
         CSK1               = 1.0
         CSK2               = CST(KV,1)
         CSK3               = CST(KV,2)
         SNK1               = 0.0
         SNK2               = SNT(KV,1)
         SNK3               = SNT(KV,2)
            DO IV              = 1,ND1
            ZVAR1R             =  WORKINP(IV,KV,  1)
            ZVAR1I             =  WORKINP(IV,KPV, 1)
            ZVAR2R             =  WORKINP(IV,KV,  2)
            ZVAR2I             =  WORKINP(IV,KPV, 2)
            ZVAR3R             =  WORKINP(IV,KV,  3)
            ZVAR3I             =  WORKINP(IV,KPV, 3)
C
            ZVAR0R1            =          +ZVAR2R      +      ZVAR1I
            ZVAR0I1            =          -ZVAR3R      +      ZVAR2I
            ZVAR0R2            =          +ZVAR1R      -  0.5*ZVAR0R1
            ZVAR0I2            =          +ZVAR3I      -  0.5*ZVAR0I1
            ZVAR0R3            = +ROOT3D2*(ZVAR2R      -      ZVAR1I)
            ZVAR0I3            = +ROOT3D2*(ZVAR3R      +      ZVAR2I)
C
            ZVAR1R1            =          +ZVAR1R      +      ZVAR0R1
            ZVAR1I1            =          +ZVAR3I      +      ZVAR0I1
            ZVAR1R2            =          +ZVAR0R2     +      ZVAR0I3
            ZVAR1I2            =          +ZVAR0I2     -      ZVAR0R3
            ZVAR1R3            =          +ZVAR0R2     -      ZVAR0I3
            ZVAR1I3            =          +ZVAR0I2     +      ZVAR0R3
C
            WORKOUT(IV,1, KV)  =    +CSK1*ZVAR1R1   -     SNK1*ZVAR1I1
            WORKOUT(IV,1, KPV) =    +CSK1*ZVAR1I1   +     SNK1*ZVAR1R1
            WORKOUT(IV,2, KV)  =    +CSK2*ZVAR1R2   -     SNK2*ZVAR1I2
            WORKOUT(IV,2, KPV) =    +CSK2*ZVAR1I2   +     SNK2*ZVAR1R2
            WORKOUT(IV,3, KV)  =    +CSK3*ZVAR1R3   -     SNK3*ZVAR1I3
            WORKOUT(IV,3, KPV) =    +CSK3*ZVAR1I3   +     SNK3*ZVAR1R3
            ENDDO
         ENDDO
      END IF
C
C
C 4.3 SPECIAL CASE FOR K = ND3/2 WITH ND3 EVEN SHOULD NOT BE NEEDED
C
      IF(N3PARITY .EQ. 0) THEN
         IERSUB   =  SUBNAM
         IERPLC   = -3
         IERROR   =  4
C
         DO IV                 = 1,ND1
         ZVAR0                 =           WORKINP(IV,ND3QP, 1)
         ZVAR1                 =           WORKINP(IV,ND3QP, 2)
         ZVAR2                 =           WORKINP(IV,ND3QP, 3)
         ZVAR3                 =         +ZVAR0    +        ZVAR1
         ZVAR3R                =         -ZVAR1    +    0.5*ZVAR0
         ZVAR3I                = +ROOT3D2*ZVAR2
C
         WORKOUT(IV,1, ND3QP)  =          +ZVAR3
         WORKOUT(IV,2, ND3QP)  =          +ZVAR3R  +        ZVAR3I
         WORKOUT(IV,3, ND3QP)  =          -ZVAR3R  +        ZVAR3I
         ENDDO
      ENDIF
C
C
C
C 5.0 SAVE THE RESULTS AND CLEAN UP
C
C 5.1 SAVE THE RESULTS
C
      IVP      = 0
      KVI      = 1
      LVI      = 1
      KVO      = 1
      LVO      = 1
      DO J     = 1,N0
      JV       = J
      DO I     = 1,M0
      IV       = I
      IVP      = IVP + 1
      IF(IVP .GT. ND1) THEN
         IVP      = 1
         LVI      = LVI + 1
         KVO      = KVO + 1
         IF(LVI .GT. ND3) THEN
            LVI       = 1
            KVI       = KVI + 1
            IF(KVI .GT. ND2) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +7
               IERROR   = KVI
               RETURN
            ENDIF
         ENDIF
         IF(KVO .GT. ND2) THEN
            KVO       = 1
            LVO       = LVO + 1
            IF(LVO .GT. ND3) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +8
               IERROR   = LVO
               RETURN
            ENDIF
         ENDIF
      ENDIF
C
      FR1    = WORKINP(IVP,LVI,KVI)
      FR2    = WORKOUT(IVP,KVO,LVO)
      IF    (SWITCH .EQ. +1) THEN
         FREAL(IV,JV)  = FR1
         FIMAG(IV,JV)  = FR2
      ELSEIF(SWITCH .EQ. -1) THEN
         FREAL(IV,JV)  = FR2
         FIMAG(IV,JV)  = FR1
      ENDIF
      ENDDO
      ENDDO
C
C
C 5.2 DEALLOCATE INTERNAL ARRAYS
C
      DEALLOCATE(WORKINP)
      DEALLOCATE(WORKOUT)
      DEALLOCATE(CST)
      DEALLOCATE(SNT)
C
C
C
C 6.0 RETURN AND END
C
      RETURN
      END

      SUBROUTINE FFTCAL4F(FREAL,FIMAG,PHASER,PHASEI,M0DIM,M0,N0DIM,N0
     &                                      ,ND1,ND2,ND3,NST,SWITCH
     &                                      ,IERSUB,IERPLC,IERROR)
C
C
C ------------------------------------------------------------------------------------
C
C CALCULATE FAST FOURIER TRANSFORM OF FREAL AND STORE IN FIMAG
C COOLEY-TOOKEY ALGORITHM FOR SPECIAL CASE WITH ND2 = 4
C
C ----------------------------------------------------------------------------------------
C
      IMPLICIT NONE
C
      CHARACTER*16 SUBNAM
      INTEGER      I,        IV,       K,        KV,       KPV
      INTEGER      IVP,      KVP,      LVP,      J,        JV
      INTEGER      KVI,      KVO,      LVI,      LVO
      INTEGER      ND3P,     ND3Q,     ND3QP
      INTEGER      NDTL0,    NDTL1,    NDTOTL0,  NDTOTL1
      INTEGER      N3PARITY
      REAL*8       ROOT2D2
      REAL*8       FR1,      FR2
      REAL*8       ZVAR0,    ZVAR1,    ZVAR2,    ZVAR3
      REAL*8       ZVAR0R,   ZVAR0I,   ZVAR1R,   ZVAR1I,
     &             ZVAR2R,   ZVAR2I
      REAL*8       ZVAR3R,   ZVAR3I,   ZVAR4R,   ZVAR4I
      REAL*8       ZVAR0R1,  ZVAR0I1,  ZVAR0R2,  ZVAR0I2
      REAL*8       ZVAR0R3,  ZVAR0I3,  ZVAR0R4,  ZVAR0I4
      REAL*8       CSK1,     CSK2,     CSK3,     CSK4
      REAL*8       SNK1,     SNK2,     SNK3,     SNK4
C
      CHARACTER*16 IERSUB
      INTEGER      M0,       N0
      INTEGER      M0DIM,    N0DIM
      INTEGER      ND1,      ND2,      ND3,      NST,      SWITCH
      INTEGER      IERPLC,   IERROR
C
      REAL*8       FREAL   (M0DIM,N0DIM),FIMAG   (M0DIM,N0DIM)
      REAL*8       PHASER  (N0DIM),      PHASEI  (N0DIM)
C
C
      REAL*8,      DIMENSION(:,:,:),ALLOCATABLE:: WORKINP
      REAL*8,      DIMENSION(:,:,:),ALLOCATABLE:: WORKOUT
      REAL*8,      DIMENSION(:,:),  ALLOCATABLE:: CST
      REAL*8,      DIMENSION(:,:),  ALLOCATABLE:: SNT
C
C
C
C 1.0 INITIALIZATION
C
      SUBNAM    = 'FFTCAL4F'
      IERSUB    = ''
      IERPLC    = 0
      IERROR    = 0
C
      ND3P      = (ND3-1)/2
      ND3Q      =  ND3/2
      ND3QP     =  ND3Q  + 1
      N3PARITY  =  ND3   - 2*(ND3 /2)
C
      IF(ND2 .NE. 4) THEN
         IERSUB   =  SUBNAM
         IERPLC   = -1
         IERROR   =  1
         RETURN
      ENDIF
C
      IF(IABS(SWITCH) .NE. +1) THEN
         IERSUB   =  SUBNAM
         IERPLC   = -2
         IERROR   =  SWITCH
         RETURN
      ENDIF
C
      ROOT2D2   = 0.707106781186547524400844362104849
C
C
C
C 2.0 STORE THE INPUT DATA IN DIFFERENT SIZE WORK ARRAYS
C
C 2.1 CHECK FOR CONSISTENT DIMENSIONS
C
      NDTOTL0   = M0*N0
      NDTOTL1   = ND1*ND2*ND3
      IF(NDTOTL0 .NE. NDTOTL1) THEN
         IERSUB   =  SUBNAM
         IERPLC   = +1
         IERROR   = NDTOTL0 - NDTOTL1
         RETURN
      ENDIF
C
      NDTL0     = ND2*ND3
      NDTL1     = N0 - NST + 1
      IF(NDTL0 .NE. NDTL1) THEN
         IERSUB   =  SUBNAM
         IERPLC   = +2
         IERROR   = NDTL0 - NDTL1
         RETURN
      ENDIF
C
      IF(NST .GT. N0) THEN
         IERSUB   =  SUBNAM
         IERPLC   = +3
         IERROR   =  NST
         RETURN
      ENDIF
C
C
C 2.2 ALLOCATE INTERNAL ARRAYS
C
      ALLOCATE (WORKINP(ND1,ND2,ND3))
      ALLOCATE (WORKOUT(ND1,ND3,ND2))
      ALLOCATE (CST    (ND3,ND2))
      ALLOCATE (SNT    (ND3,ND2))
C
C
C 2.3 FILL IN THE INTERNAL ARRAYS
C
C 2.3.1 FILL IN THE WORK ARRAYS WITH DIMENSIONS (ND1,ND2,ND3) AND (ND1,ND3,ND2)
C
      IVP      = 0
      KVI      = 1
      LVI      = 1
      KVO      = 1
      LVO      = 1
      DO J     = 1,N0
      JV       = J
      DO I     = 1,M0
      IV       = I
      IVP      = IVP + 1
      IF(IVP .GT. ND1) THEN
         IVP      = 1
         KVI      = KVI + 1
         LVO      = LVO + 1
         IF(KVI .GT. ND2) THEN
            KVI       = 1
            LVI       = LVI + 1
            IF(LVI .GT. ND3) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +4
               IERROR   = LVI
               RETURN
            ENDIF
         ENDIF
         IF(LVO .GT. ND3) THEN
            LVO       = 1
            KVO       = KVO + 1
            IF(KVO .GT. ND2) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +5
               IERROR   = KVO
               RETURN
            ENDIF
         ENDIF
      ENDIF
C
      IF    (SWITCH .EQ. +1) THEN
         FR1    = FREAL(IV,JV)
      ELSEIF(SWITCH .EQ. -1) THEN
         FR1    = FIMAG(IV,JV)
      ENDIF
      WORKINP(IVP,KVI,LVI)  = FR1
      WORKOUT(IVP,LVO,KVO)  = 0.0
      ENDDO
      ENDDO
C
C 2.3.2 FILL IN THE PHASE ARRAYS FROM K = NST TO N0
C
      LVP      = 0
      KVP      = 1
      DO J     = NST,N0
      JV       = J
      LVP      = LVP + 1
      IF(LVP .GT. ND3) THEN
         LVP      = 1
         KVP      = KVP + 1
         IF(KVP .GT. ND2) THEN
            IERSUB   =  SUBNAM
            IERPLC   = +6
            IERROR   = KVP
            RETURN
         ENDIF
      ENDIF
      CST(LVP,KVP)  = PHASER(JV)
      SNT(LVP,KVP)  = PHASEI(JV)
      ENDDO
C
C
C
C 3.0 KV = 1: SET FIRST ELEMENT OF MIDDLE INDEX
C
      DO IV              = 1,ND1
      ZVAR0R             =           WORKINP(IV,1,1)  +
     &                               WORKINP(IV,3,1)
      ZVAR0I             =           WORKINP(IV,2,1)  +
     &                               WORKINP(IV,4,1)
      ZVAR1R             =           WORKINP(IV,1,1)  -
     &                               WORKINP(IV,3,1)
      ZVAR1I             =           WORKINP(IV,2,1)  -
     &                               WORKINP(IV,4,1)
C
      WORKOUT(IV,1,1)    =          +ZVAR0R           +   ZVAR0I
      WORKOUT(IV,1,2)    =          +ZVAR1R
      WORKOUT(IV,1,3)    =          +ZVAR0R           -   ZVAR0I
      WORKOUT(IV,1,4)    =          -ZVAR1I
      ENDDO
C
C
C
C 4.0 FOR KV = 2, ND3P + 1
C
C 4.1 FOR ND1 > (ND3-1)/2
C
      IF    (ND1 .GT. ND3P) THEN
         DO IV              = 1,ND1
            DO K               = 1,ND3P
            KV                 = K   + 1
            KPV                = ND3 - K + 1
            CSK1               = 1.0
            CSK2               = CST(KV,1)
            CSK3               = CST(KV,2)
            CSK4               = CST(KV,3)
            SNK1               = 0.0
            SNK2               = SNT(KV,1)
            SNK3               = SNT(KV,2)
            SNK4               = SNT(KV,3)
C
            ZVAR1R             = CSK1*WORKINP(IV,1, KV)  -
     &                                          SNK1*WORKINP(IV,1, KPV) 
            ZVAR1I             = CSK1*WORKINP(IV,1, KPV) +
     &                                          SNK1*WORKINP(IV,1, KV)
            ZVAR2R             = CSK2*WORKINP(IV,2, KV)  -
     &                                          SNK2*WORKINP(IV,2, KPV)
            ZVAR2I             = CSK2*WORKINP(IV,2, KPV) +
     &                                          SNK2*WORKINP(IV,2, KV)
            ZVAR3R             = CSK3*WORKINP(IV,3, KV)  -
     &                                          SNK3*WORKINP(IV,3, KPV)
            ZVAR3I             = CSK3*WORKINP(IV,3, KPV) +
     &                                          SNK3*WORKINP(IV,3, KV)
            ZVAR4R             = CSK4*WORKINP(IV,4, KV)  -
     &                                          SNK4*WORKINP(IV,4, KPV)
            ZVAR4I             = CSK4*WORKINP(IV,4, KPV) +
     &                                          SNK4*WORKINP(IV,4, KV)
C
            ZVAR0R1            =      ZVAR1R   +      ZVAR3R
            ZVAR0I1            =      ZVAR1I   +      ZVAR3I
            ZVAR0R2            =      ZVAR2R   +      ZVAR4R
            ZVAR0I2            =      ZVAR2I   +      ZVAR4I
            ZVAR0R3            =      ZVAR1R   -      ZVAR3R
            ZVAR0I3            =      ZVAR1I   -      ZVAR3I
            ZVAR0R4            =      ZVAR2R   -      ZVAR4R
            ZVAR0I4            =      ZVAR2I   -      ZVAR4I
C
            WORKOUT(IV,KV, 1)  =      ZVAR0R1  +      ZVAR0R2
            WORKOUT(IV,KPV,1)  =      ZVAR0R3  -      ZVAR0I4
            WORKOUT(IV,KV, 2)  =      ZVAR0R3  +      ZVAR0I4
            WORKOUT(IV,KPV,2)  =      ZVAR0R1  -      ZVAR0R2
            WORKOUT(IV,KV, 3)  =     -ZVAR0I1  +      ZVAR0I2
            WORKOUT(IV,KPV,3)  =      ZVAR0I3  -      ZVAR0R4
            WORKOUT(IV,KV, 4)  =     -ZVAR0I3  -      ZVAR0R4
            WORKOUT(IV,KPV,4)  =      ZVAR0I1  +      ZVAR0I2
            ENDDO
         ENDDO
C
C
C 4.2 FOR ND1 > (ND3-1)/2
C
      ELSEIF(ND1 .LE. ND3P) THEN
         DO K               = 1,ND3P
         KV                 = K   + 1
         KPV                = ND3 - K + 1
         CSK1               = 1.0
         CSK2               = CST(KV,1)
         CSK3               = CST(KV,2)
         CSK4               = CST(KV,3)
         SNK1               = 0.0
         SNK2               = SNT(KV,1)
         SNK3               = SNT(KV,2)
         SNK4               = SNT(KV,3)
            DO IV              = 1,ND1
            ZVAR1R             = CSK1*WORKINP(IV,1, KV)  -
     &                                          SNK1*WORKINP(IV,1,KPV)
            ZVAR1I             = CSK1*WORKINP(IV,1, KPV) +
     &                                          SNK1*WORKINP(IV,1,KV)
            ZVAR2R             = CSK2*WORKINP(IV,2, KV)  -
     &                                          SNK2*WORKINP(IV,2,KPV)
            ZVAR2I             = CSK2*WORKINP(IV,2, KPV) +
     &                                          SNK2*WORKINP(IV,2,KV)
            ZVAR3R             = CSK3*WORKINP(IV,3, KV)  -
     &                                          SNK3*WORKINP(IV,3,KPV)
            ZVAR3I             = CSK3*WORKINP(IV,3, KPV) +
     &                                          SNK3*WORKINP(IV,3,KV)
            ZVAR4R             = CSK4*WORKINP(IV,4, KV)  -
     &                                          SNK4*WORKINP(IV,4,KPV)
            ZVAR4I             = CSK4*WORKINP(IV,4, KPV) +
     &                                          SNK4*WORKINP(IV,4,KV)
C
            ZVAR0R1            =      ZVAR1R   +      ZVAR3R
            ZVAR0I1            =      ZVAR1I   +      ZVAR3I
            ZVAR0R2            =      ZVAR2R   +      ZVAR4R
            ZVAR0I2            =      ZVAR2I   +      ZVAR4I
            ZVAR0R3            =      ZVAR1R   -      ZVAR3R
            ZVAR0I3            =      ZVAR1I   -      ZVAR3I
            ZVAR0R4            =      ZVAR2R   -      ZVAR4R
            ZVAR0I4            =      ZVAR2I   -      ZVAR4I
C
            WORKOUT(IV,KV, 1)  =      ZVAR0R1  +      ZVAR0R2
            WORKOUT(IV,KPV,1)  =      ZVAR0R3  -      ZVAR0I4
            WORKOUT(IV,KV, 2)  =      ZVAR0R3  +      ZVAR0I4
            WORKOUT(IV,KPV,2)  =      ZVAR0R1  -      ZVAR0R2
            WORKOUT(IV,KV, 3)  =     -ZVAR0I1  +      ZVAR0I2
            WORKOUT(IV,KPV,3)  =      ZVAR0I3  -      ZVAR0R4
            WORKOUT(IV,KV, 4)  =     -ZVAR0I3  -      ZVAR0R4
            WORKOUT(IV,KPV,4)  =      ZVAR0I1  +      ZVAR0I2
            ENDDO
         ENDDO
      END IF
C
C
C 4.3 SPECIAL CASE FOR K = ND3/2 WITH ND3 EVEN
C
      IF(N3PARITY .EQ. 0) THEN
         DO IV                 = 1,ND1
         ZVAR0                 =           WORKINP(IV,1,ND3QP)
         ZVAR1                 =  ROOT2D2*(WORKINP(IV,2,ND3QP) -
     &                                     WORKINP(IV,4,ND3QP))
         ZVAR2                 =  ROOT2D2*(WORKINP(IV,2,ND3QP) +
     &                                     WORKINP(IV,4,ND3QP))
         ZVAR3                 =           WORKINP(IV,3,ND3QP)
C
         WORKOUT(IV,ND3QP, 1)  =          +ZVAR0               + ZVAR1
         WORKOUT(IV,ND3QP, 2)  =          +ZVAR0               - ZVAR1
         WORKOUT(IV,ND3QP, 3)  =          -ZVAR2               + ZVAR3
         WORKOUT(IV,ND3QP, 4)  =          -ZVAR2               - ZVAR3
         ENDDO
      ENDIF
C
C
C
C 5.0 SAVE THE RESULTS AND CLEAN UP
C
C 5.1 SAVE THE RESULTS
C
      IVP      = 0
      KVI      = 1
      LVI      = 1
      KVO      = 1
      LVO      = 1
      DO J     = 1,N0
      JV       = J
      DO I     = 1,M0
      IV       = I
      IVP      = IVP + 1
      IF(IVP .GT. ND1) THEN
         IVP      = 1
         KVI      = KVI + 1
         LVO      = LVO + 1
         IF(KVI .GT. ND2) THEN
            KVI       = 1
            LVI       = LVI + 1
            IF(LVI .GT. ND3) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +7
               IERROR   = LVI
               RETURN
            ENDIF
         ENDIF
         IF(LVO .GT. ND3) THEN
            LVO       = 1
            KVO       = KVO + 1
            IF(KVO .GT. ND2) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +8
               IERROR   = KVO
               RETURN
            ENDIF
         ENDIF
      ENDIF
C
      FR1    = WORKINP(IVP,KVI,LVI)
      FR2    = WORKOUT(IVP,LVO,KVO)
      IF    (SWITCH .EQ. +1) THEN
         FREAL(IV,JV)  = FR1
         FIMAG(IV,JV)  = FR2
      ELSEIF(SWITCH .EQ. -1) THEN
         FREAL(IV,JV)  = FR2
         FIMAG(IV,JV)  = FR1
      ENDIF
      ENDDO
      ENDDO
C
C
C 5.2 DEALLOCATE INTERNAL ARRAYS
C
      DEALLOCATE(WORKINP)
      DEALLOCATE(WORKOUT)
      DEALLOCATE(CST)
      DEALLOCATE(SNT)
C
C
C
C 6.0 RETURN AND END
C
      RETURN
      END
      SUBROUTINE FFTCAL4I(FREAL,FIMAG,PHASER,PHASEI,M0DIM,M0,N0DIM,N0
     &                                      ,ND1,ND2,ND3,NST,SWITCH
     &                                      ,IERSUB,IERPLC,IERROR)
C
C
C ------------------------------------------------------------------------------------
C
C CALCULATE INVERSE FAST FOURIER TRANSFORM OF FREAL AND STORE IN FIMAG
C COOLEY-TOOKEY ALGORITHM FOR SPECIAL CASE WITH ND2 = 4
C
C ----------------------------------------------------------------------------------------
C
      IMPLICIT NONE
C
      CHARACTER*16 SUBNAM
      INTEGER      I,        IV,       K,        KV,       KPV
      INTEGER      IVP,      KVP,      LVP,      J,        JV
      INTEGER      KVI,      KVO,      LVI,      LVO
      INTEGER      ND3P,     ND3Q,     ND3QP
      INTEGER      NDTL0,    NDTL1,    NDTOTL0,  NDTOTL1
      INTEGER      N3PARITY
      REAL*8       ROOT2D2
      REAL*8       FR1,      FR2
      REAL*8       ZVAR0,    ZVAR1,    ZVAR2,    ZVAR3
      REAL*8       ZVAR0R,   ZVAR0I,   ZVAR1R,   ZVAR1I,
     &             ZVAR2R,   ZVAR2I
      REAL*8       ZVAR3R,   ZVAR3I,   ZVAR4R,   ZVAR4I
      REAL*8       ZVAR0R1,  ZVAR0I1,  ZVAR0R2,  ZVAR0I2
      REAL*8       ZVAR0R3,  ZVAR0I3,  ZVAR0R4,  ZVAR0I4
      REAL*8       CSK1,     CSK2,     CSK3,     CSK4
      REAL*8       SNK1,     SNK2,     SNK3,     SNK4
C
      CHARACTER*16 IERSUB
      INTEGER      M0,       N0
      INTEGER      M0DIM,    N0DIM
      INTEGER      ND1,      ND2,      ND3,      NST,      SWITCH
      INTEGER      IERPLC,   IERROR
C
      REAL*8       FREAL   (M0DIM,N0DIM),FIMAG   (M0DIM,N0DIM)
      REAL*8       PHASER  (N0DIM),      PHASEI  (N0DIM)
C
C
      REAL*8,      DIMENSION(:,:,:),ALLOCATABLE:: WORKINP
      REAL*8,      DIMENSION(:,:,:),ALLOCATABLE:: WORKOUT
      REAL*8,      DIMENSION(:,:),  ALLOCATABLE:: CST
      REAL*8,      DIMENSION(:,:),  ALLOCATABLE:: SNT
C
C
C
C 1.0 INITIALIZATION
C
      SUBNAM    = 'FFTCAL4I'
      IERSUB    = ''
      IERPLC    = 0
      IERROR    = 0
C
      ND3P      = (ND3-1)/2
      ND3Q      =  ND3/2
      ND3QP     =  ND3Q  + 1
      N3PARITY  =  ND3   - 2*(ND3 /2)
C
      IF(ND2 .NE. 4) THEN
         IERSUB   =  SUBNAM
         IERPLC   = -1
         IERROR   =  1
         RETURN
      ENDIF
C
      IF(IABS(SWITCH) .NE. +1) THEN
         IERSUB   =  SUBNAM
         IERPLC   = -2
         IERROR   =  SWITCH
         RETURN
      ENDIF
C
      ROOT2D2   = 0.707106781186547524400844362104849
C
C
C
C 2.0 STORE THE INPUT DATA IN DIFFERENT SIZE WORK ARRAYS
C
C 2.1 CHECK FOR CONSISTENT DIMENSIONS
C
      NDTOTL0   = M0*N0
      NDTOTL1   = ND1*ND2*ND3
      IF(NDTOTL0 .NE. NDTOTL1) THEN
         IERSUB   =  SUBNAM
         IERPLC   = +1
         IERROR   = NDTOTL0 - NDTOTL1
         RETURN
      ENDIF
C
      NDTL0     = ND2*ND3
      NDTL1     = N0 - NST + 1
      IF(NDTL0 .NE. NDTL1) THEN
         IERSUB   =  SUBNAM
         IERPLC   = +2
         IERROR   = NDTL0 - NDTL1
         RETURN
      ENDIF
C
      IF(NST .GT. N0) THEN
         IERSUB   =  SUBNAM
         IERPLC   = +3
         IERROR   =  NST
         RETURN
      ENDIF
C
C
C 2.2 ALLOCATE INTERNAL ARRAYS
C
      ALLOCATE (WORKINP(ND1,ND3,ND2))
      ALLOCATE (WORKOUT(ND1,ND2,ND3))
      ALLOCATE (CST    (ND3,ND2))
      ALLOCATE (SNT    (ND3,ND2))
C
C
C 2.3 FILL IN THE INTERNAL ARRAYS
C
C 2.3.1 FILL IN THE WORK ARRAYS WITH DIMENSIONS (ND1,ND3,ND2) AND (ND1,ND2,ND3)
C
      IVP      = 0
      KVI      = 1
      LVI      = 1
      KVO      = 1
      LVO      = 1
      DO J     = 1,N0
      JV       = J
      DO I     = 1,M0
      IV       = I
      IVP      = IVP + 1
      IF(IVP .GT. ND1) THEN
         IVP      = 1
         LVI      = LVI + 1
         KVO      = KVO + 1
         IF(LVI .GT. ND3) THEN
            LVI       = 1
            KVI       = KVI + 1
            IF(KVI .GT. ND2) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +4
               IERROR   = KVI
               RETURN
            ENDIF
         ENDIF
         IF(KVO .GT. ND2) THEN
            KVO       = 1
            LVO       = LVO + 1
            IF(LVO .GT. ND3) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +5
               IERROR   = LVO
               RETURN
            ENDIF
         ENDIF
      ENDIF
C
      IF    (SWITCH .EQ. +1) THEN
         FR1    = FREAL(IV,JV)
      ELSEIF(SWITCH .EQ. -1) THEN
         FR1    = FIMAG(IV,JV)
      ENDIF
      WORKINP(IVP,LVI,KVI)  = FR1
      WORKOUT(IVP,KVO,LVO)  = 0.0
      ENDDO
      ENDDO
C
C 2.3.2 FILL IN THE PHASE ARRAYS FROM K = NST TO N0
C
      LVP      = 0
      KVP      = 1
      DO J     = NST,N0
      JV       = J
      LVP      = LVP + 1
      IF(LVP .GT. ND3) THEN
         LVP      = 1
         KVP      = KVP + 1
         IF(KVP .GT. ND2) THEN
            IERSUB   =  SUBNAM
            IERPLC   = +6
            IERROR   = KVP
            RETURN
         ENDIF
      ENDIF
      CST(LVP,KVP)  = PHASER(JV)
      SNT(LVP,KVP)  = PHASEI(JV)
      ENDDO
C
C
C
C 3.0 KV = 1: SET FIRST ELEMENT OF MIDDLE INDEX
C
      DO IV              = 1,ND1
      ZVAR0R             =           WORKINP(IV,1,1)  +
     &                               WORKINP(IV,1,3)
      ZVAR0I             =           WORKINP(IV,1,1)  -
     &                               WORKINP(IV,1,3)
      ZVAR1R             =           WORKINP(IV,1,2)
      ZVAR1I             =           WORKINP(IV,1,4)
C
      WORKOUT(IV,1,1)    =          +ZVAR0R           +   ZVAR1R
      WORKOUT(IV,2,1)    =          +ZVAR0I           +   ZVAR1I
      WORKOUT(IV,3,1)    =          +ZVAR0R           -   ZVAR1R
      WORKOUT(IV,4,1)    =          +ZVAR0I           -   ZVAR1I
      ENDDO
C
C
C
C 4.0 FOR KV = 2, ND3P + 1
C
C 4.1 FOR ND1 > (ND3-1)/2
C
      IF    (ND1 .GT. ND3P) THEN
         DO IV              = 1,ND1
            DO K               = 1,ND3P
            KV                 = K   + 1
            KPV                = ND3 - K + 1
            CSK1               = 1.0
            CSK2               = CST(KV,1)
            CSK3               = CST(KV,2)
            CSK4               = CST(KV,3)
            SNK1               = 0.0
            SNK2               = SNT(KV,1)
            SNK3               = SNT(KV,2)
            SNK4               = SNT(KV,3)
C
            ZVAR1R             = +WORKINP(IV,KV,  1)   +
     &                            WORKINP(IV,KPV, 2)
            ZVAR1I             = -WORKINP(IV,KV,  3)   +
     &                            WORKINP(IV,KPV, 4)
            ZVAR2R             = +WORKINP(IV,KV,  2)   +
     &                            WORKINP(IV,KPV, 1)
            ZVAR2I             = -WORKINP(IV,KV,  4)   +
     &                            WORKINP(IV,KPV, 3)
            ZVAR3R             = +WORKINP(IV,KV,  1)   -
     &                            WORKINP(IV,KPV, 2)
            ZVAR3I             = +WORKINP(IV,KV,  3)   +
     &                            WORKINP(IV,KPV, 4)
            ZVAR4R             = +WORKINP(IV,KV,  2)   -
     &                            WORKINP(IV,KPV, 1)
            ZVAR4I             = +WORKINP(IV,KV,  4)   +
     &                            WORKINP(IV,KPV, 3)
C
            ZVAR0R1            =      ZVAR1R   +      ZVAR2R
            ZVAR0I1            =      ZVAR1I   +      ZVAR2I
            ZVAR0R2            =      ZVAR3R   +      ZVAR4I
            ZVAR0I2            =      ZVAR3I   -      ZVAR4R
            ZVAR0R3            =      ZVAR1R   -      ZVAR2R
            ZVAR0I3            =      ZVAR1I   -      ZVAR2I
            ZVAR0R4            =      ZVAR3R   -      ZVAR4I
            ZVAR0I4            =      ZVAR3I   +      ZVAR4R
C
            WORKOUT(IV,1, KV)  =    +CSK1*ZVAR0R1   -     SNK1*ZVAR0I1
            WORKOUT(IV,1, KPV) =    +CSK1*ZVAR0I1   +     SNK1*ZVAR0R1
            WORKOUT(IV,2, KV)  =    +CSK2*ZVAR0R2   -     SNK2*ZVAR0I2
            WORKOUT(IV,2, KPV) =    +CSK2*ZVAR0I2   +     SNK2*ZVAR0R2
            WORKOUT(IV,3, KV)  =    +CSK3*ZVAR0R3   -     SNK3*ZVAR0I3
            WORKOUT(IV,3, KPV) =    +CSK3*ZVAR0I3   +     SNK3*ZVAR0R3
            WORKOUT(IV,4, KV)  =    +CSK4*ZVAR0R4   -     SNK4*ZVAR0I4
            WORKOUT(IV,4, KPV) =    +CSK4*ZVAR0I4   +     SNK4*ZVAR0R4
            ENDDO
         ENDDO
C
C
C 4.2 FOR ND1 > (ND3-1)/2
C
      ELSEIF(ND1 .LE. ND3P) THEN
         DO K               = 1,ND3P
         KV                 = K   + 1
         KPV                = ND3 - K + 1
         CSK1               = 1.0
         CSK2               = CST(KV,1)
         CSK3               = CST(KV,2)
         CSK4               = CST(KV,3)
         SNK1               = 0.0
         SNK2               = SNT(KV,1)
         SNK3               = SNT(KV,2)
         SNK4               = SNT(KV,3)
            DO IV              = 1,ND1
            ZVAR1R             = +WORKINP(IV,KV,  1)   +
     &                            WORKINP(IV,KPV, 2)
            ZVAR1I             = -WORKINP(IV,KV,  3)   +
     &                            WORKINP(IV,KPV, 4)
            ZVAR2R             = +WORKINP(IV,KV,  2)   +
     &                            WORKINP(IV,KPV, 1)
            ZVAR2I             = -WORKINP(IV,KV,  4)   +
     &                            WORKINP(IV,KPV, 3)
            ZVAR3R             = +WORKINP(IV,KV,  1)   -
     &                            WORKINP(IV,KPV, 2)
            ZVAR3I             = +WORKINP(IV,KV,  3)   +
     &                            WORKINP(IV,KPV, 4)
            ZVAR4R             = +WORKINP(IV,KV,  2)   -
     &                            WORKINP(IV,KPV, 1)
            ZVAR4I             = +WORKINP(IV,KV,  4)   +
     &                            WORKINP(IV,KPV, 3)
C
            ZVAR0R1            =      ZVAR1R   +      ZVAR2R
            ZVAR0I1            =      ZVAR1I   +      ZVAR2I
            ZVAR0R2            =      ZVAR3R   +      ZVAR4I
            ZVAR0I2            =      ZVAR3I   -      ZVAR4R
            ZVAR0R3            =      ZVAR1R   -      ZVAR2R
            ZVAR0I3            =      ZVAR1I   -      ZVAR2I
            ZVAR0R4            =      ZVAR3R   -      ZVAR4I
            ZVAR0I4            =      ZVAR3I   +      ZVAR4R
C
            WORKOUT(IV,1, KV)  =    +CSK1*ZVAR0R1   -     SNK1*ZVAR0I1
            WORKOUT(IV,1, KPV) =    +CSK1*ZVAR0I1   +     SNK1*ZVAR0R1
            WORKOUT(IV,2, KV)  =    +CSK2*ZVAR0R2   -     SNK2*ZVAR0I2
            WORKOUT(IV,2, KPV) =    +CSK2*ZVAR0I2   +     SNK2*ZVAR0R2
            WORKOUT(IV,3, KV)  =    +CSK3*ZVAR0R3   -     SNK3*ZVAR0I3
            WORKOUT(IV,3, KPV) =    +CSK3*ZVAR0I3   +     SNK3*ZVAR0R3
            WORKOUT(IV,4, KV)  =    +CSK4*ZVAR0R4   -     SNK4*ZVAR0I4
            WORKOUT(IV,4, KPV) =    +CSK4*ZVAR0I4   +     SNK4*ZVAR0R4
            ENDDO
         ENDDO
      END IF
C
C
C 4.3 SPECIAL CASE FOR K = ND3/2 WITH ND3 EVEN
C
      IF(N3PARITY .EQ. 0) THEN
         DO IV                 = 1,ND1
         ZVAR0                 =  WORKINP(IV,ND3QP, 1)
         ZVAR1                 =  WORKINP(IV,ND3QP, 2)
         ZVAR2                 =  WORKINP(IV,ND3QP, 3)
         ZVAR3                 =  WORKINP(IV,ND3QP, 4)
         ZVAR0R                =          +ZVAR0           + ZVAR1
         ZVAR0I                =          -ZVAR2           + ZVAR3
         ZVAR1R                =          +ZVAR0           - ZVAR1
         ZVAR1I                =          +ZVAR2           + ZVAR3
C
         WORKOUT(IV,1, ND3QP)  =          +ZVAR0R
         WORKOUT(IV,2, ND3QP)  = +ROOT2D2*(ZVAR1R          + ZVAR1I)
         WORKOUT(IV,3, ND3QP)  =          +ZVAR0I
         WORKOUT(IV,4, ND3QP)  = -ROOT2D2*(ZVAR1R          - ZVAR1I)
         ENDDO
      ENDIF
C
C
C
C 5.0 SAVE THE RESULTS AND CLEAN UP
C
C 5.1 SAVE THE RESULTS
C
      IVP      = 0
      KVI      = 1
      LVI      = 1
      KVO      = 1
      LVO      = 1
      DO J     = 1,N0
      JV       = J
      DO I     = 1,M0
      IV       = I
      IVP      = IVP + 1
      IF(IVP .GT. ND1) THEN
         IVP      = 1
         LVI      = LVI + 1
         KVO      = KVO + 1
         IF(LVI .GT. ND3) THEN
            LVI       = 1
            KVI       = KVI + 1
            IF(KVI .GT. ND2) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +7
               IERROR   = KVI
               RETURN
            ENDIF
         ENDIF
         IF(KVO .GT. ND2) THEN
            KVO       = 1
            LVO       = LVO + 1
            IF(LVO .GT. ND3) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +8
               IERROR   = LVO
               RETURN
            ENDIF
         ENDIF
      ENDIF
C
      FR1    = WORKINP(IVP,LVI,KVI)
      FR2    = WORKOUT(IVP,KVO,LVO)
      IF    (SWITCH .EQ. +1) THEN
         FREAL(IV,JV)  = FR1
         FIMAG(IV,JV)  = FR2
      ELSEIF(SWITCH .EQ. -1) THEN
         FREAL(IV,JV)  = FR2
         FIMAG(IV,JV)  = FR1
      ENDIF
      ENDDO
      ENDDO
C
C
C 5.2 DEALLOCATE INTERNAL ARRAYS
C
      DEALLOCATE(WORKINP)
      DEALLOCATE(WORKOUT)
      DEALLOCATE(CST)
      DEALLOCATE(SNT)
C
C
C
C 6.0 RETURN AND END
C
      RETURN
      END

      SUBROUTINE FFTCAL5F(FREAL,FIMAG,PHASER,PHASEI,M0DIM,M0,N0DIM,N0
     &                                      ,ND1,ND2,ND3,NST,SWITCH
     &                                      ,IERSUB,IERPLC,IERROR)
C
C
C ------------------------------------------------------------------------------------
C
C CALCULATE FAST FOURIER TRANSFORM OF FREAL AND STORE IN FIMAG
C COOLEY-TOOKEY ALGORITHM FOR SPECIAL CASE WITH ND2 = 5
C
C ----------------------------------------------------------------------------------------
C
      IMPLICIT NONE
C
      CHARACTER*16 SUBNAM
      INTEGER      I,        IV,       K,        KV,       KPV
      INTEGER      IVP,      KVP,      LVP,      J,        JV
      INTEGER      KVI,      KVO,      LVI,      LVO
      INTEGER      ND3P,     ND3Q,     ND3QP
      INTEGER      NDTL0,    NDTL1,    NDTOTL0,  NDTOTL1
      INTEGER      N3PARITY
      REAL*8       SIN36,    SIN72,    SINRAT,   ROOT5D4
      REAL*8       FR1,      FR2
      REAL*8       ZVAR0,    ZVAR3,    ZVAR4
      REAL*8       ZVAR0R,   ZVAR0I,   ZVAR1R,   ZVAR1I,
     &             ZVAR2R,   ZVAR2I
      REAL*8       ZVAR3R,   ZVAR3I,   ZVAR4R,   ZVAR4I
      REAL*8       ZVAR0R1,  ZVAR0I1,  ZVAR0R2,  ZVAR0I2
      REAL*8       ZVAR1R1,  ZVAR1I1,  ZVAR1R2,  ZVAR1I2
      REAL*8       ZVAR2R0,  ZVAR2I0,  ZVAR2R1,  ZVAR2I1,
     &             ZVAR2R2,  ZVAR2I2
      REAL*8       ZVAR3R1,  ZVAR3I1,  ZVAR3R2,  ZVAR3I2
      REAL*8       ZVAR4R1,  ZVAR4I1,  ZVAR4R2,  ZVAR4I2
      REAL*8       CSK1,     CSK2,     CSK3,     CSK4,     CSK5
      REAL*8       SNK1,     SNK2,     SNK3,     SNK4,     SNK5
C
      CHARACTER*16 IERSUB
      INTEGER      M0,       N0
      INTEGER      M0DIM,    N0DIM
      INTEGER      ND1,      ND2,      ND3,      NST,      SWITCH
      INTEGER      IERPLC,   IERROR
C
      REAL*8       FREAL   (M0DIM,N0DIM),FIMAG   (M0DIM,N0DIM)
      REAL*8       PHASER  (N0DIM),      PHASEI  (N0DIM)
C
C
      REAL*8,      DIMENSION(:,:,:),ALLOCATABLE:: WORKINP
      REAL*8,      DIMENSION(:,:,:),ALLOCATABLE:: WORKOUT
      REAL*8,      DIMENSION(:,:),  ALLOCATABLE:: CST
      REAL*8,      DIMENSION(:,:),  ALLOCATABLE:: SNT
C
C
C
C 1.0 INITIALIZATION
C
      SUBNAM    = 'FFTCAL5F'
      IERSUB    = ''
      IERPLC    = 0
      IERROR    = 0
C
      ND3P      = (ND3-1)/2
      ND3Q      =  ND3/2
      ND3QP     =  ND3Q  + 1
      N3PARITY  =  ND3   - 2*(ND3 /2)
C
      IF(ND2 .NE. 5) THEN
         IERSUB   =  SUBNAM
         IERPLC   = -1
         IERROR   =  1
         RETURN
      ENDIF
C
      IF(IABS(SWITCH) .NE. +1) THEN
         IERSUB   =  SUBNAM
         IERPLC   = -2
         IERROR   =  SWITCH
         RETURN
      ENDIF
C
      ROOT5D4   = 0.559016994374947424102293417182819
      SIN72     = 0.951056516295153572116439333379382
      SIN36     = 0.587785252292473137103456792829093
      SINRAT    = SIN36/SIN72
C
C
C
C 2.0 STORE THE INPUT DATA IN DIFFERENT SIZE WORK ARRAYS
C
C 2.1 CHECK FOR CONSISTENT DIMENSIONS
C
      NDTOTL0   = M0*N0
      NDTOTL1   = ND1*ND2*ND3
      IF(NDTOTL0 .NE. NDTOTL1) THEN
         IERSUB   =  SUBNAM
         IERPLC   = +1
         IERROR   = NDTOTL0 - NDTOTL1
         RETURN
      ENDIF
C
      NDTL0     = ND2*ND3
      NDTL1     = N0 - NST + 1
      IF(NDTL0 .NE. NDTL1) THEN
         IERSUB   =  SUBNAM
         IERPLC   = +2
         IERROR   = NDTL0 - NDTL1
         RETURN
      ENDIF
C
      IF(NST .GT. N0) THEN
         IERSUB   =  SUBNAM
         IERPLC   = +3
         IERROR   =  NST
         RETURN
      ENDIF
C
C
C 2.2 ALLOCATE INTERNAL ARRAYS
C
      ALLOCATE (WORKINP(ND1,ND2,ND3))
      ALLOCATE (WORKOUT(ND1,ND3,ND2))
      ALLOCATE (CST    (ND3,ND2))
      ALLOCATE (SNT    (ND3,ND2))
C
C
C 2.3 FILL IN THE INTERNAL ARRAYS
C
C 2.3.1 FILL IN THE WORK ARRAYS WITH DIMENSIONS (ND1,ND2,ND3) AND (ND1,ND3,ND2)
C
      IVP      = 0
      KVI      = 1
      LVI      = 1
      KVO      = 1
      LVO      = 1
      DO J     = 1,N0
      JV       = J
      DO I     = 1,M0
      IV       = I
      IVP      = IVP + 1
      IF(IVP .GT. ND1) THEN
         IVP      = 1
         KVI      = KVI + 1
         LVO      = LVO + 1
         IF(KVI .GT. ND2) THEN
            KVI       = 1
            LVI       = LVI + 1
            IF(LVI .GT. ND3) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +4
               IERROR   = LVI
               RETURN
            ENDIF
         ENDIF
         IF(LVO .GT. ND3) THEN
            LVO       = 1
            KVO       = KVO + 1
            IF(KVO .GT. ND2) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +5
               IERROR   = KVO
               RETURN
            ENDIF
         ENDIF
      ENDIF
C
      IF    (SWITCH .EQ. +1) THEN
         FR1    = FREAL(IV,JV)
      ELSEIF(SWITCH .EQ. -1) THEN
         FR1    = FIMAG(IV,JV)
      ENDIF
      WORKINP(IVP,KVI,LVI)  = FR1
      WORKOUT(IVP,LVO,KVO)  = 0.0
      ENDDO
      ENDDO
C
C 2.3.2 FILL IN THE PHASE ARRAYS FROM K = NST TO N0
C
      LVP      = 0
      KVP      = 1
      DO J     = NST,N0
      JV       = J
      LVP      = LVP + 1
      IF(LVP .GT. ND3) THEN
         LVP      = 1
         KVP      = KVP + 1
         IF(KVP .GT. ND2) THEN
            IERSUB   =  SUBNAM
            IERPLC   = +6
            IERROR   = KVP
            RETURN
         ENDIF
      ENDIF
      CST(LVP,KVP)  = PHASER(JV)
      SNT(LVP,KVP)  = PHASEI(JV)
      ENDDO
C
C
C
C 3.0 KV = 1: SET FIRST ELEMENT OF MIDDLE INDEX
C
      DO IV              = 1,ND1
      ZVAR0              =           WORKINP(IV,1,1)
      ZVAR0R             =           WORKINP(IV,2,1)  +
     &                               WORKINP(IV,5,1)
      ZVAR1R             =           WORKINP(IV,3,1)  +
     &                               WORKINP(IV,4,1)
      ZVAR0I             =  SIN72  *(WORKINP(IV,2,1)  -
     &                               WORKINP(IV,5,1))
      ZVAR1I             =  SIN72  *(WORKINP(IV,3,1)  -
     &                               WORKINP(IV,4,1))
      ZVAR2R             =           ZVAR0R           +      ZVAR1R
      ZVAR2I             =  ROOT5D4*(ZVAR0R           -      ZVAR1R)
      ZVAR3              =           ZVAR0            - 0.25*ZVAR2R
C
      WORKOUT(IV,1,1)    =           ZVAR0            +      ZVAR2R
      WORKOUT(IV,1,2)    =           ZVAR3            +      ZVAR2I
      WORKOUT(IV,1,3)    =           ZVAR3            -      ZVAR2I
      WORKOUT(IV,1,4)    = -SINRAT  *ZVAR0I           +      ZVAR1I
      WORKOUT(IV,1,5)    = -SINRAT  *ZVAR1I           -      ZVAR0I
      ENDDO
C
C
C
C 4.0 FOR KV = 2, ND3P + 1
C
C 4.1 FOR ND1 > (ND3-1)/2
C
      IF    (ND1 .GT. ND3P) THEN
         DO IV              = 1,ND1
            DO K               = 1,ND3P
            KV                 = K   + 1
            KPV                = ND3 - K + 1
            CSK1               = 1.0
            CSK2               = CST(KV,1)
            CSK3               = CST(KV,2)
            CSK4               = CST(KV,3)
            CSK5               = CST(KV,4)
            SNK1               = 0.0
            SNK2               = SNT(KV,1)
            SNK3               = SNT(KV,2)
            SNK4               = SNT(KV,3)
            SNK5               = SNT(KV,4)
C
            ZVAR0R             = CSK1*WORKINP(IV,1, KV)  -
     &                                          SNK1*WORKINP(IV,1, KPV)
            ZVAR0I             = CSK1*WORKINP(IV,1, KPV) +
     &                                          SNK1*WORKINP(IV,1, KV)
            ZVAR1R             = CSK2*WORKINP(IV,2, KV)  -
     &                                          SNK2*WORKINP(IV,2, KPV)
            ZVAR1I             = CSK2*WORKINP(IV,2, KPV) +
     &                                          SNK2*WORKINP(IV,2, KV)
            ZVAR2R             = CSK3*WORKINP(IV,3, KV)  -
     &                                          SNK3*WORKINP(IV,3, KPV)
            ZVAR2I             = CSK3*WORKINP(IV,3, KPV) +
     &                                          SNK3*WORKINP(IV,3, KV)
            ZVAR3R             = CSK4*WORKINP(IV,4, KV)  -
     &                                          SNK4*WORKINP(IV,4, KPV)
            ZVAR3I             = CSK4*WORKINP(IV,4, KPV) +
     &                                          SNK4*WORKINP(IV,4, KV)
            ZVAR4R             = CSK5*WORKINP(IV,5, KV)  -
     &                                          SNK5*WORKINP(IV,5, KPV)
            ZVAR4I             = CSK5*WORKINP(IV,5, KPV) +
     &                                          SNK5*WORKINP(IV,5, KV)
C
            ZVAR0R1            =           ZVAR1R   +      ZVAR4R
            ZVAR0I1            =           ZVAR1I   +      ZVAR4I
            ZVAR0R2            =           ZVAR2R   +      ZVAR3R
            ZVAR0I2            =           ZVAR2I   +      ZVAR3I
C
            ZVAR1R1            =  SIN72  *(ZVAR1R   -      ZVAR4R)
            ZVAR1I1            =  SIN72  *(ZVAR1I   -      ZVAR4I)
            ZVAR1R2            =  SIN72  *(ZVAR2R   -      ZVAR3R)
            ZVAR1I2            =  SIN72  *(ZVAR2I   -      ZVAR3I)
C
            ZVAR2R0            =           ZVAR0R1  +      ZVAR0R2
            ZVAR2I0            =           ZVAR0I1  +      ZVAR0I2
            ZVAR2R1            =  ROOT5D4*(ZVAR0R1  -      ZVAR0R2)
            ZVAR2I1            =  ROOT5D4*(ZVAR0I1  -      ZVAR0I2)
            ZVAR2R2            =           ZVAR0R   - 0.25*ZVAR2R0
            ZVAR2I2            =           ZVAR0I   - 0.25*ZVAR2I0
C
            ZVAR3R1            =           ZVAR2R2  +      ZVAR2R1
            ZVAR3I1            =           ZVAR2I2  +      ZVAR2I1
            ZVAR3R2            =           ZVAR2R2  -      ZVAR2R1
            ZVAR3I2            =           ZVAR2I2  -      ZVAR2I1
C
            ZVAR4R1            =  SINRAT  *ZVAR1R2  +      ZVAR1R1
            ZVAR4I1            =  SINRAT  *ZVAR1I2  +      ZVAR1I1
            ZVAR4R2            =  SINRAT  *ZVAR1R1  -      ZVAR1R2
            ZVAR4I2            =  SINRAT  *ZVAR1I1  -      ZVAR1I2
C
            WORKOUT(IV,KV, 1)  =           ZVAR0R   +      ZVAR2R0
            WORKOUT(IV,KPV,1)  =           ZVAR3R1  -      ZVAR4I1
            WORKOUT(IV,KV, 2)  =           ZVAR3R1  +      ZVAR4I1
            WORKOUT(IV,KPV,2)  =           ZVAR3R2  -      ZVAR4I2
            WORKOUT(IV,KV, 3)  =           ZVAR3R2  +      ZVAR4I2
            WORKOUT(IV,KPV,3)  =           ZVAR3I2  -      ZVAR4R2
            WORKOUT(IV,KV, 4)  =          -ZVAR3I2  -      ZVAR4R2
            WORKOUT(IV,KPV,4)  =           ZVAR3I1  -      ZVAR4R1
            WORKOUT(IV,KV, 5)  =          -ZVAR3I1  -      ZVAR4R1
            WORKOUT(IV,KPV,5)  =           ZVAR0I   +      ZVAR2I0
            ENDDO
         ENDDO
C
C
C 4.2 FOR ND1 > (ND3-1)/2
C
      ELSEIF(ND1 .LE. ND3P) THEN
         DO K               = 1,ND3P
         KV                 = K   + 1
         KPV                = ND3 - K + 1
         CSK1               = 1.0
         CSK2               = CST(KV,1)
         CSK3               = CST(KV,2)
         CSK4               = CST(KV,3)
         CSK5               = CST(KV,4)
         SNK1               = 0.0
         SNK2               = SNT(KV,1)
         SNK3               = SNT(KV,2)
         SNK4               = SNT(KV,3)
         SNK5               = SNT(KV,4)
            DO IV              = 1,ND1
            ZVAR0R             = CSK1*WORKINP(IV,1, KV)  -
     &                                          SNK1*WORKINP(IV,1, KPV)
            ZVAR0I             = CSK1*WORKINP(IV,1, KPV) +
     &                                          SNK1*WORKINP(IV,1, KV)
            ZVAR1R             = CSK2*WORKINP(IV,2, KV)  -
     &                                          SNK2*WORKINP(IV,2, KPV)
            ZVAR1I             = CSK2*WORKINP(IV,2, KPV) +
     &                                          SNK2*WORKINP(IV,2, KV)
            ZVAR2R             = CSK3*WORKINP(IV,3, KV)  -
     &                                          SNK3*WORKINP(IV,3, KPV)
            ZVAR2I             = CSK3*WORKINP(IV,3, KPV) +
     &                                          SNK3*WORKINP(IV,3, KV)
            ZVAR3R             = CSK4*WORKINP(IV,4, KV)  -
     &                                          SNK4*WORKINP(IV,4, KPV)
            ZVAR3I             = CSK4*WORKINP(IV,4, KPV) +
     &                                          SNK4*WORKINP(IV,4, KV)
            ZVAR4R             = CSK5*WORKINP(IV,5, KV)  -
     &                                          SNK5*WORKINP(IV,5, KPV)
            ZVAR4I             = CSK5*WORKINP(IV,5, KPV) +
     &                                          SNK5*WORKINP(IV,5, KV)
C
            ZVAR0R1            =           ZVAR1R   +      ZVAR4R
            ZVAR0I1            =           ZVAR1I   +      ZVAR4I
            ZVAR0R2            =           ZVAR2R   +      ZVAR3R
            ZVAR0I2            =           ZVAR2I   +      ZVAR3I
C
            ZVAR1R1            =  SIN72  *(ZVAR1R   -      ZVAR4R)
            ZVAR1I1            =  SIN72  *(ZVAR1I   -      ZVAR4I)
            ZVAR1R2            =  SIN72  *(ZVAR2R   -      ZVAR3R)
            ZVAR1I2            =  SIN72  *(ZVAR2I   -      ZVAR3I)
C
            ZVAR2R0            =           ZVAR0R1  +      ZVAR0R2
            ZVAR2I0            =           ZVAR0I1  +      ZVAR0I2
            ZVAR2R1            =  ROOT5D4*(ZVAR0R1  -      ZVAR0R2)
            ZVAR2I1            =  ROOT5D4*(ZVAR0I1  -      ZVAR0I2)
            ZVAR2R2            =           ZVAR0R   - 0.25*ZVAR2R0
            ZVAR2I2            =           ZVAR0I   - 0.25*ZVAR2I0
C
            ZVAR3R1            =           ZVAR2R2  +      ZVAR2R1
            ZVAR3I1            =           ZVAR2I2  +      ZVAR2I1
            ZVAR3R2            =           ZVAR2R2  -      ZVAR2R1
            ZVAR3I2            =           ZVAR2I2  -      ZVAR2I1
C
            ZVAR4R1            =  SINRAT  *ZVAR1R2  +      ZVAR1R1
            ZVAR4I1            =  SINRAT  *ZVAR1I2  +      ZVAR1I1
            ZVAR4R2            =  SINRAT  *ZVAR1R1  -      ZVAR1R2
            ZVAR4I2            =  SINRAT  *ZVAR1I1  -      ZVAR1I2
C
            WORKOUT(IV,KV, 1)  =           ZVAR0R   +      ZVAR2R0
            WORKOUT(IV,KPV,1)  =           ZVAR3R1  -      ZVAR4I1
            WORKOUT(IV,KV, 2)  =           ZVAR3R1  +      ZVAR4I1
            WORKOUT(IV,KPV,2)  =           ZVAR3R2  -      ZVAR4I2
            WORKOUT(IV,KV, 3)  =           ZVAR3R2  +      ZVAR4I2
            WORKOUT(IV,KPV,3)  =           ZVAR3I2  -      ZVAR4R2
            WORKOUT(IV,KV, 4)  =          -ZVAR3I2  -      ZVAR4R2
            WORKOUT(IV,KPV,4)  =           ZVAR3I1  -      ZVAR4R1
            WORKOUT(IV,KV, 5)  =          -ZVAR3I1  -      ZVAR4R1
            WORKOUT(IV,KPV,5)  =           ZVAR0I   +      ZVAR2I0
            ENDDO
         ENDDO
      END IF
C
C
C 4.3 SPECIAL CASE FOR K = ND3/2 WITH ND3 EVEN SHOULD NOT BE NEEDED
C
      IF(N3PARITY .EQ. 0) THEN
         IERSUB   =  SUBNAM
         IERPLC   = -3
         IERROR   =  4
C
         DO IV                = 1,ND1
         ZVAR0                =           WORKINP(IV,1, ND3QP)
         ZVAR1R               =           WORKINP(IV,2, ND3QP) -
     &                                    WORKINP(IV,5, ND3QP)
         ZVAR1I               = +SIN72*  (WORKINP(IV,2, ND3QP) +
     &                                    WORKINP(IV,5, ND3QP))
         ZVAR2R               =           WORKINP(IV,3, ND3QP) -
     &                                    WORKINP(IV,4, ND3QP)
         ZVAR2I               = +SIN72*  (WORKINP(IV,3, ND3QP) +
     &                                    WORKINP(IV,4, ND3QP))
         ZVAR3R               =          +ZVAR2R     -         ZVAR1R
         ZVAR3I               = +ROOT5D4*(ZVAR2R     +         ZVAR1R)
         ZVAR4                =          +ZVAR0      -    0.25*ZVAR3R 
C
         WORKOUT(IV,ND3QP, 1) =          +ZVAR4      +         ZVAR3I
         WORKOUT(IV,ND3QP, 2) =          +ZVAR4      -         ZVAR3I
         WORKOUT(IV,ND3QP, 3) =          +ZVAR0      +         ZVAR3R
         WORKOUT(IV,ND3QP, 4) =          -ZVAR1I     +  SINRAT*ZVAR2I
         WORKOUT(IV,ND3QP, 5) = -SINRAT*  ZVAR1I     -         ZVAR2I
         ENDDO
      ENDIF
C
C
C
C 5.0 SAVE THE RESULTS AND CLEAN UP
C
C 5.1 SAVE THE RESULTS
C
      IVP      = 0
      KVI      = 1
      LVI      = 1
      KVO      = 1
      LVO      = 1
      DO J     = 1,N0
      JV       = J
      DO I     = 1,M0
      IV       = I
      IVP      = IVP + 1
      IF(IVP .GT. ND1) THEN
         IVP      = 1
         KVI      = KVI + 1
         LVO      = LVO + 1
         IF(KVI .GT. ND2) THEN
            KVI       = 1
            LVI       = LVI + 1
            IF(LVI .GT. ND3) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +7
               IERROR   = LVI
               RETURN
            ENDIF
         ENDIF
         IF(LVO .GT. ND3) THEN
            LVO       = 1
            KVO       = KVO + 1
            IF(KVO .GT. ND2) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +8
               IERROR   = KVO
               RETURN
            ENDIF
         ENDIF
      ENDIF
C
      FR1    = WORKINP(IVP,KVI,LVI)
      FR2    = WORKOUT(IVP,LVO,KVO)
      IF    (SWITCH .EQ. +1) THEN
         FREAL(IV,JV)  = FR1
         FIMAG(IV,JV)  = FR2
      ELSEIF(SWITCH .EQ. -1) THEN
         FREAL(IV,JV)  = FR2
         FIMAG(IV,JV)  = FR1
      ENDIF
      ENDDO
      ENDDO
C
C
C 5.2 DEALLOCATE INTERNAL ARRAYS
C
      DEALLOCATE(WORKINP)
      DEALLOCATE(WORKOUT)
      DEALLOCATE(CST)
      DEALLOCATE(SNT)
C
C
C
C 6.0 RETURN AND END
C
      RETURN
      END
      SUBROUTINE FFTCAL5I(FREAL,FIMAG,PHASER,PHASEI,M0DIM,M0,N0DIM,N0
     &                                      ,ND1,ND2,ND3,NST,SWITCH
     &                                      ,IERSUB,IERPLC,IERROR)
C
C
C ------------------------------------------------------------------------------------
C
C CALCULATE INVERSE FAST FOURIER TRANSFORM OF FREAL AND STORE IN FIMAG
C COOLEY-TOOKEY ALGORITHM FOR SPECIAL CASE WITH ND2 = 5
C
C ----------------------------------------------------------------------------------------
C
      IMPLICIT NONE
C
      CHARACTER*16 SUBNAM
      INTEGER      I,        IV,       K,        KV,       KPV
      INTEGER      IVP,      KVP,      LVP,      J,        JV
      INTEGER      KVI,      KVO,      LVI,      LVO
      INTEGER      ND3P,     ND3Q,     ND3QP
      INTEGER      NDTL0,    NDTL1,    NDTOTL0,  NDTOTL1
      INTEGER      N3PARITY
      REAL*8       SIN36,    SIN72,    SINRAT,   ROOT5D4
      REAL*8       FR1,      FR2
      REAL*8       ZVAR0,    ZVAR1,    ZVAR2,    ZVAR3,    ZVAR4,
     &             ZVAR5
      REAL*8       ZVAR0R,   ZVAR0I,   ZVAR1R,   ZVAR1I,
     &             ZVAR2R,   ZVAR2I
      REAL*8       ZVAR3R,   ZVAR3I,   ZVAR4R,   ZVAR4I,
     &             ZVAR5R,   ZVAR5I
      REAL*8       ZVAR0R1,  ZVAR0I1,  ZVAR0R2,  ZVAR0I2,
     &             ZVAR0R3,  ZVAR0I3,  ZVAR0R4,  ZVAR0I4
      REAL*8       ZVAR1R1,  ZVAR1I1,  ZVAR1R2,  ZVAR1I2,
     &             ZVAR1R3,  ZVAR1I3
      REAL*8       ZVAR2R0,  ZVAR2I0,  ZVAR2R1,  ZVAR2I1,
     &             ZVAR2R2,  ZVAR2I2,  ZVAR2R3,  ZVAR2I3,
     &             ZVAR2R4,  ZVAR2I4
      REAL*8       CSK1,     CSK2,     CSK3,     CSK4,     CSK5
      REAL*8       SNK1,     SNK2,     SNK3,     SNK4,     SNK5
C
      CHARACTER*16 IERSUB
      INTEGER      M0,       N0
      INTEGER      M0DIM,    N0DIM
      INTEGER      ND1,      ND2,      ND3,      NST,      SWITCH
      INTEGER      IERPLC,   IERROR
C
      REAL*8       FREAL   (M0DIM,N0DIM),FIMAG   (M0DIM,N0DIM)
      REAL*8       PHASER  (N0DIM),      PHASEI  (N0DIM)
C
C
      REAL*8,      DIMENSION(:,:,:),ALLOCATABLE:: WORKINP
      REAL*8,      DIMENSION(:,:,:),ALLOCATABLE:: WORKOUT
      REAL*8,      DIMENSION(:,:),  ALLOCATABLE:: CST
      REAL*8,      DIMENSION(:,:),  ALLOCATABLE:: SNT
C
C
C
C 1.0 INITIALIZATION
C
      SUBNAM    = 'FFTCAL5I'
      IERSUB    = ''
      IERPLC    = 0
      IERROR    = 0
C
      ND3P      = (ND3-1)/2
      ND3Q      =  ND3/2
      ND3QP     =  ND3Q + 1
      N3PARITY  =  ND3  - 2*(ND3 /2)
C
      IF(ND2 .NE. 5) THEN
         IERSUB   =  SUBNAM
         IERPLC   = -1
         IERROR   =  1
         RETURN
      ENDIF
C
      IF(IABS(SWITCH) .NE. +1) THEN
         IERSUB   =  SUBNAM
         IERPLC   = -2
         IERROR   =  SWITCH
         RETURN
      ENDIF
C
      ROOT5D4   = 0.559016994374947424102293417182819
      SIN72     = 0.951056516295153572116439333379382
      SIN36     = 0.587785252292473137103456792829093
      SINRAT    = SIN36/SIN72
C
C
C
C 2.0 STORE THE INPUT DATA IN DIFFERENT SIZE WORK ARRAYS
C
C 2.1 CHECK FOR CONSISTENT DIMENSIONS
C
      NDTOTL0   = M0*N0
      NDTOTL1   = ND1*ND2*ND3
      IF(NDTOTL0 .NE. NDTOTL1) THEN
         IERSUB   =  SUBNAM
         IERPLC   = +1
         IERROR   = NDTOTL0 - NDTOTL1
         RETURN
      ENDIF
C
      NDTL0     = ND2*ND3
      NDTL1     = N0 - NST + 1
      IF(NDTL0 .NE. NDTL1) THEN
         IERSUB   =  SUBNAM
         IERPLC   = +2
         IERROR   = NDTL0 - NDTL1
         RETURN
      ENDIF
C
      IF(NST .GT. N0) THEN
         IERSUB   =  SUBNAM
         IERPLC   = +3
         IERROR   =  NST
         RETURN
      ENDIF
C
C
C 2.2 ALLOCATE INTERNAL ARRAYS
C
      ALLOCATE (WORKINP(ND1,ND3,ND2))
      ALLOCATE (WORKOUT(ND1,ND2,ND3))
      ALLOCATE (CST    (ND3,ND2))
      ALLOCATE (SNT    (ND3,ND2))
C
C
C 2.3 FILL IN THE INTERNAL ARRAYS
C
C 2.3.1 FILL IN THE WORK ARRAYS WITH DIMENSIONS (ND1,ND3,ND2) AND (ND1,ND2,ND3)
C
      IVP      = 0
      KVI      = 1
      LVI      = 1
      KVO      = 1
      LVO      = 1
      DO J     = 1,N0
      JV       = J
      DO I     = 1,M0
      IV       = I
      IVP      = IVP + 1
      IF(IVP .GT. ND1) THEN
         IVP      = 1
         LVI      = LVI + 1
         KVO      = KVO + 1
         IF(LVI .GT. ND3) THEN
            LVI       = 1
            KVI       = KVI + 1
            IF(KVI .GT. ND2) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +4
               IERROR   = KVI
               RETURN
            ENDIF
         ENDIF
         IF(KVO .GT. ND2) THEN
            KVO       = 1
            LVO       = LVO + 1
            IF(LVO .GT. ND3) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +5
               IERROR   = LVO
               RETURN
            ENDIF
         ENDIF
      ENDIF
C
      IF    (SWITCH .EQ. +1) THEN
         FR1    = FREAL(IV,JV)
      ELSEIF(SWITCH .EQ. -1) THEN
         FR1    = FIMAG(IV,JV)
      ENDIF
      WORKINP(IVP,LVI,KVI)  = FR1
      WORKOUT(IVP,KVO,LVO)  = 0.0
      ENDDO
      ENDDO
C
C 2.3.2 FILL IN THE PHASE ARRAYS FROM K = NST TO N0
C
      LVP      = 0
      KVP      = 1
      DO J     = NST,N0
      JV       = J
      LVP      = LVP + 1
      IF(LVP .GT. ND3) THEN
         LVP      = 1
         KVP      = KVP + 1
         IF(KVP .GT. ND2) THEN
            IERSUB   =  SUBNAM
            IERPLC   = +6
            IERROR   = KVP
            RETURN
         ENDIF
      ENDIF
      CST(LVP,KVP)  = PHASER(JV)
      SNT(LVP,KVP)  = PHASEI(JV)
      ENDDO
C
C
C
C 3.0 KV = 1: SET FIRST ELEMENT OF MIDDLE INDEX
C
      DO IV              = 1,ND1
      ZVAR0              =           WORKINP(IV,1,  1)
      ZVAR1              =           WORKINP(IV,1,  2)
      ZVAR2              =           WORKINP(IV,1,  3)
      ZVAR3              =  SIN72*   WORKINP(IV,1,  5)
      ZVAR4              =  SIN72*   WORKINP(IV,1,  4)
      ZVAR0R             =           ZVAR1     + ZVAR2
      ZVAR0I             =  ROOT5D4*(ZVAR1     - ZVAR2)
      ZVAR5              =           ZVAR0     - 0.25*ZVAR0R
      ZVAR1R             =          +ZVAR5     +      ZVAR0I
      ZVAR1I             =          +ZVAR5     -      ZVAR0I
      ZVAR2R             =  SINRAT*  ZVAR4     +      ZVAR3
      ZVAR2I             =  SINRAT*  ZVAR3     -      ZVAR4
C
      WORKOUT(IV,1, 1)   = +ZVAR0   + ZVAR0R
      WORKOUT(IV,2, 1)   = +ZVAR1R  + ZVAR2R
      WORKOUT(IV,3, 1)   = +ZVAR1I  + ZVAR2I
      WORKOUT(IV,4, 1)   = +ZVAR1I  - ZVAR2I
      WORKOUT(IV,5, 1)   = +ZVAR1R  - ZVAR2R
      ENDDO
C
C
C
C 4.0 FOR KV = 2, ND3P + 1
C
C 4.1 FOR ND1 > (ND3-1)/2
C
      IF    (ND1 .GT. ND3P) THEN
         DO IV              = 1,ND1
            DO K               = 1,ND3P
            KV                 = K   + 1
            KPV                = ND3 - K + 1
            CSK1               = 1.0
            CSK2               = CST(KV,1)
            CSK3               = CST(KV,2)
            CSK4               = CST(KV,3)
            CSK5               = CST(KV,4)
            SNK1               = 0.0
            SNK2               = SNT(KV,1)
            SNK3               = SNT(KV,2)
            SNK4               = SNT(KV,3)
            SNK5               = SNT(KV,4)
C
            ZVAR1R             = WORKINP(IV,KV,  1)
            ZVAR1I             = WORKINP(IV,KPV, 1)
            ZVAR2R             = WORKINP(IV,KV,  2)
            ZVAR2I             = WORKINP(IV,KPV, 2)
            ZVAR3R             = WORKINP(IV,KV,  3)
            ZVAR3I             = WORKINP(IV,KPV, 3)
            ZVAR4R             = WORKINP(IV,KV,  4)
            ZVAR4I             = WORKINP(IV,KPV, 4)
            ZVAR5R             = WORKINP(IV,KV,  5)
            ZVAR5I             = WORKINP(IV,KPV, 5)
C
            ZVAR0R1            =          +ZVAR2R   +          ZVAR1I
            ZVAR0I1            =  SIN72  *(ZVAR2R   -          ZVAR1I)
            ZVAR0R2            =          -ZVAR5R   +          ZVAR4I
            ZVAR0I2            =  SIN72  *(ZVAR5R   +          ZVAR4I)
            ZVAR0R3            =          +ZVAR3R   +          ZVAR2I
            ZVAR0I3            =  SIN72  *(ZVAR3R   -          ZVAR2I)
            ZVAR0R4            =          -ZVAR4R   +          ZVAR3I
            ZVAR0I4            =  SIN72  *(ZVAR4R   +          ZVAR3I)
C
            ZVAR1R1            =          +ZVAR0R1  +          ZVAR0R3
            ZVAR1I1            =          +ZVAR0R2  +          ZVAR0R4
            ZVAR1R2            =  ROOT5D4*(ZVAR0R1  -          ZVAR0R3)
            ZVAR1I2            =  ROOT5D4*(ZVAR0R2  -          ZVAR0R4)
            ZVAR1R3            =          +ZVAR1R   -     0.25*ZVAR1R1
            ZVAR1I3            =          +ZVAR5I   -     0.25*ZVAR1I1
C
            ZVAR2R0            =           ZVAR1R
            ZVAR2I0            =           ZVAR5I
            ZVAR2R1            =           ZVAR1R3  +          ZVAR1R2
            ZVAR2I1            =           ZVAR1I3  +          ZVAR1I2
            ZVAR2R2            =           ZVAR1R3  -          ZVAR1R2
            ZVAR2I2            =           ZVAR1I3  -          ZVAR1I2
            ZVAR2R3            =           ZVAR0I1  +   SINRAT*ZVAR0I3
            ZVAR2I3            =           ZVAR0I2  +   SINRAT*ZVAR0I4
            ZVAR2R4            =  SINRAT  *ZVAR0I1  -          ZVAR0I3
            ZVAR2I4            =  SINRAT  *ZVAR0I2  -          ZVAR0I4
C
            ZVAR1R             =           ZVAR2R0  +          ZVAR1R1
            ZVAR1I             =           ZVAR2I0  +          ZVAR1I1
            ZVAR2R             =           ZVAR2R1  +          ZVAR2I3
            ZVAR2I             =           ZVAR2I1  -          ZVAR2R3
            ZVAR3R             =           ZVAR2R2  +          ZVAR2I4
            ZVAR3I             =           ZVAR2I2  -          ZVAR2R4
            ZVAR4R             =           ZVAR2R2  -          ZVAR2I4
            ZVAR4I             =           ZVAR2I2  +          ZVAR2R4
            ZVAR5R             =           ZVAR2R1  -          ZVAR2I3
            ZVAR5I             =           ZVAR2I1  +          ZVAR2R3
C
            WORKOUT(IV,1, KV)  =    +CSK1*ZVAR1R    -     SNK1*ZVAR1I
            WORKOUT(IV,1, KPV) =    +CSK1*ZVAR1I    +     SNK1*ZVAR1R
            WORKOUT(IV,2, KV)  =    +CSK2*ZVAR2R    -     SNK2*ZVAR2I
            WORKOUT(IV,2, KPV) =    +CSK2*ZVAR2I    +     SNK2*ZVAR2R
            WORKOUT(IV,3, KV)  =    +CSK3*ZVAR3R    -     SNK3*ZVAR3I
            WORKOUT(IV,3, KPV) =    +CSK3*ZVAR3I    +     SNK3*ZVAR3R
            WORKOUT(IV,4, KV)  =    +CSK4*ZVAR4R    -     SNK4*ZVAR4I
            WORKOUT(IV,4, KPV) =    +CSK4*ZVAR4I    +     SNK4*ZVAR4R
            WORKOUT(IV,5, KV)  =    +CSK5*ZVAR5R    -     SNK5*ZVAR5I
            WORKOUT(IV,5, KPV) =    +CSK5*ZVAR5I    +     SNK5*ZVAR5R
            ENDDO
         ENDDO
C
C
C 4.2 FOR ND1 > (ND3-1)/2
C
      ELSEIF(ND1 .LE. ND3P) THEN
         DO K               = 1,ND3P
         KV                 = K   + 1
         KPV                = ND3 - K + 1
         CSK1               = 1.0
         CSK2               = CST(KV,1)
         CSK3               = CST(KV,2)
         CSK4               = CST(KV,3)
         CSK5               = CST(KV,4)
         SNK1               = 0.0
         SNK2               = SNT(KV,1)
         SNK3               = SNT(KV,2)
         SNK4               = SNT(KV,3)
         SNK5               = SNT(KV,4)
            DO IV              = 1,ND1
            ZVAR1R             = WORKINP(IV,KV,  1)
            ZVAR1I             = WORKINP(IV,KPV, 1)
            ZVAR2R             = WORKINP(IV,KV,  2)
            ZVAR2I             = WORKINP(IV,KPV, 2)
            ZVAR3R             = WORKINP(IV,KV,  3)
            ZVAR3I             = WORKINP(IV,KPV, 3)
            ZVAR4R             = WORKINP(IV,KV,  4)
            ZVAR4I             = WORKINP(IV,KPV, 4)
            ZVAR5R             = WORKINP(IV,KV,  5)
            ZVAR5I             = WORKINP(IV,KPV, 5)
C
            ZVAR0R1            =          +ZVAR2R   +          ZVAR1I
            ZVAR0I1            =  SIN72  *(ZVAR2R   -          ZVAR1I)
            ZVAR0R2            =          -ZVAR5R   +          ZVAR4I
            ZVAR0I2            =  SIN72  *(ZVAR5R   +          ZVAR4I)
            ZVAR0R3            =          +ZVAR3R   +          ZVAR2I
            ZVAR0I3            =  SIN72  *(ZVAR3R   -          ZVAR2I)
            ZVAR0R4            =          -ZVAR4R   +          ZVAR3I
            ZVAR0I4            =  SIN72  *(ZVAR4R   +          ZVAR3I)
C
            ZVAR1R1            =          +ZVAR0R1  +          ZVAR0R3
            ZVAR1I1            =          +ZVAR0R2  +          ZVAR0R4
            ZVAR1R2            =  ROOT5D4*(ZVAR0R1  -          ZVAR0R3)
            ZVAR1I2            =  ROOT5D4*(ZVAR0R2  -          ZVAR0R4)
            ZVAR1R3            =          +ZVAR1R   -     0.25*ZVAR1R1
            ZVAR1I3            =          +ZVAR5I   -     0.25*ZVAR1I1
C
            ZVAR2R0            =           ZVAR1R
            ZVAR2I0            =           ZVAR5I
            ZVAR2R1            =           ZVAR1R3  +          ZVAR1R2
            ZVAR2I1            =           ZVAR1I3  +          ZVAR1I2
            ZVAR2R2            =           ZVAR1R3  -          ZVAR1R2
            ZVAR2I2            =           ZVAR1I3  -          ZVAR1I2
            ZVAR2R3            =           ZVAR0I1  +   SINRAT*ZVAR0I3
            ZVAR2I3            =           ZVAR0I2  +   SINRAT*ZVAR0I4
            ZVAR2R4            =  SINRAT  *ZVAR0I1  -          ZVAR0I3
            ZVAR2I4            =  SINRAT  *ZVAR0I2  -          ZVAR0I4
C
            ZVAR1R             =           ZVAR2R0  +          ZVAR1R1
            ZVAR1I             =           ZVAR2I0  +          ZVAR1I1
            ZVAR2R             =           ZVAR2R1  +          ZVAR2I3
            ZVAR2I             =           ZVAR2I1  -          ZVAR2R3
            ZVAR3R             =           ZVAR2R2  +          ZVAR2I4
            ZVAR3I             =           ZVAR2I2  -          ZVAR2R4
            ZVAR4R             =           ZVAR2R2  -          ZVAR2I4
            ZVAR4I             =           ZVAR2I2  +          ZVAR2R4
            ZVAR5R             =           ZVAR2R1  -          ZVAR2I3
            ZVAR5I             =           ZVAR2I1  +          ZVAR2R3
C
            WORKOUT(IV,1, KV)  =    +CSK1*ZVAR1R    -     SNK1*ZVAR1I
            WORKOUT(IV,1, KPV) =    +CSK1*ZVAR1I    +     SNK1*ZVAR1R
            WORKOUT(IV,2, KV)  =    +CSK2*ZVAR2R    -     SNK2*ZVAR2I
            WORKOUT(IV,2, KPV) =    +CSK2*ZVAR2I    +     SNK2*ZVAR2R
            WORKOUT(IV,3, KV)  =    +CSK3*ZVAR3R    -     SNK3*ZVAR3I
            WORKOUT(IV,3, KPV) =    +CSK3*ZVAR3I    +     SNK3*ZVAR3R
            WORKOUT(IV,4, KV)  =    +CSK4*ZVAR4R    -     SNK4*ZVAR4I
            WORKOUT(IV,4, KPV) =    +CSK4*ZVAR4I    +     SNK4*ZVAR4R
            WORKOUT(IV,5, KV)  =    +CSK5*ZVAR5R    -     SNK5*ZVAR5I
            WORKOUT(IV,5, KPV) =    +CSK5*ZVAR5I    +     SNK5*ZVAR5R
            ENDDO
         ENDDO
      ENDIF
C
C
C 4.3 SPECIAL CASE FOR K = ND3/2 WITH ND3 EVEN
C
      IF(N3PARITY .EQ. 0) THEN
         IERSUB   =  SUBNAM
         IERPLC   = -3
         IERROR   =  4
C
         DO IV                = 1,ND1
         ZVAR1                =           WORKINP(IV,ND3QP,3 )
         ZVAR1R               =           WORKINP(IV,ND3QP,1 ) +
     &                                    WORKINP(IV,ND3QP,2 )
         ZVAR1I               = +ROOT5D4*(WORKINP(IV,ND3QP,1 ) -
     &                                    WORKINP(IV,ND3QP,2 ))
         ZVAR2                =          -ZVAR1                +
     &                                                     0.25*ZVAR1R
         ZVAR3R               = +SIN36*   WORKINP(IV,ND3QP,5 ) +
     &                           SIN72*   WORKINP(IV,ND3QP,4 )
         ZVAR3I               = -SIN36*   WORKINP(IV,ND3QP,4 ) +
     &                           SIN72*   WORKINP(IV,ND3QP,5 )
         ZVAR4R               =          +ZVAR2      +         ZVAR1I
         ZVAR4I               =          +ZVAR2      -         ZVAR1I
C
         WORKOUT(IV,1,ND3QP)  =          +ZVAR1R     +         ZVAR1
         WORKOUT(IV,2,ND3QP)  =          +ZVAR3R     +         ZVAR4R
         WORKOUT(IV,3,ND3QP)  =          +ZVAR3I     -         ZVAR4I
         WORKOUT(IV,4,ND3QP)  =          +ZVAR3I     +         ZVAR4I
         WORKOUT(IV,5,ND3QP)  =          +ZVAR3R     -         ZVAR4R
         ENDDO
      ENDIF
C
C
C
C 5.0 SAVE THE RESULTS AND CLEAN UP
C
C 5.1 SAVE THE RESULTS
C
      IVP      = 0
      KVI      = 1
      LVI      = 1
      KVO      = 1
      LVO      = 1
      DO J     = 1,N0
      JV       = J
      DO I     = 1,M0
      IV       = I
      IVP      = IVP + 1
      IF(IVP .GT. ND1) THEN
         IVP      = 1
         LVI      = LVI + 1
         KVO      = KVO + 1
         IF(LVI .GT. ND3) THEN
            LVI       = 1
            KVI       = KVI + 1
            IF(KVI .GT. ND2) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +7
               IERROR   = KVI
               RETURN
            ENDIF
         ENDIF
         IF(KVO .GT. ND2) THEN
            KVO       = 1
            LVO       = LVO + 1
            IF(LVO .GT. ND3) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +8
               IERROR   = LVO
               RETURN
            ENDIF
         ENDIF
      ENDIF
C
      FR1    = WORKINP(IVP,LVI,KVI)
      FR2    = WORKOUT(IVP,KVO,LVO)
      IF    (SWITCH .EQ. +1) THEN
         FREAL(IV,JV)  = FR1
         FIMAG(IV,JV)  = FR2
      ELSEIF(SWITCH .EQ. -1) THEN
         FREAL(IV,JV)  = FR2
         FIMAG(IV,JV)  = FR1
      ENDIF
      ENDDO
      ENDDO
C
C
C 5.2 DEALLOCATE INTERNAL ARRAYS
C
      DEALLOCATE(WORKINP)
      DEALLOCATE(WORKOUT)
      DEALLOCATE(CST)
      DEALLOCATE(SNT)
C
C
C
C 6.0 RETURN AND END
C
      RETURN
      END
      SUBROUTINE FFTCAL6F(FREAL,FIMAG,PHASER,PHASEI,M0DIM,M0,N0DIM,N0
     &                                      ,ND1,ND2,ND3,NST,SWITCH
     &                                      ,IERSUB,IERPLC,IERROR)
C
C
C ------------------------------------------------------------------------------------
C
C CALCULATE FAST FOURIER TRANSFORM OF FREAL AND STORE IN FIMAG
C COOLEY-TOOKEY ALGORITHM FOR SPECIAL CASE WITH ND2 = 6
C
C ----------------------------------------------------------------------------------------
C
      IMPLICIT NONE
C
      CHARACTER*16 IERSUB
      INTEGER      I,        IV,       K,        KV,       KPV
      INTEGER      IVP,      KVP,      LVP,      J,        JV
      INTEGER      KVI,      KVO,      LVI,      LVO
      INTEGER      ND3P,     ND3Q,     ND3QP
      INTEGER      NDTL0,    NDTL1,    NDTOTL0,  NDTOTL1
      INTEGER      N3PARITY
      REAL*8       ROOT3D2
      REAL*8       FR1,      FR2
      REAL*8       ZVAR0,    ZVAR1,    ZVAR2,    ZVAR3
      REAL*8       ZVAR0R,   ZVAR0I,   ZVAR1R,   ZVAR1I,   ZVAR2R,
     &             ZVAR2I,   ZVAR3R,   ZVAR3I,   ZVAR4R,   ZVAR4I,
     &             ZVAR5R,   ZVAR5I
      REAL*8       ZVAR0R1,  ZVAR0I1,  ZVAR0R2,  ZVAR0I2,
     &             ZVAR0R3,  ZVAR0I3
      REAL*8       ZVAR1R1,  ZVAR1I1,  ZVAR1R2,  ZVAR1I2,
     &             ZVAR1R3,  ZVAR1I3
      REAL*8       ZVAR2R1,  ZVAR2I1,  ZVAR2R2,  ZVAR2I2,
     &             ZVAR2R3,  ZVAR2I3
      REAL*8       ZVAR3R1,  ZVAR3I1,  ZVAR3R2,  ZVAR3I2,
     &             ZVAR3R3,  ZVAR3I3
C
      REAL*8       CSK1,     CSK2,     CSK3,     CSK4,
     &             CSK5,     CSK6
      REAL*8       SNK1,     SNK2,     SNK3,     SNK4,
     &             SNK5,     SNK6
C
      CHARACTER*16 SUBNAM
      INTEGER      M0,       N0
      INTEGER      M0DIM,    N0DIM
      INTEGER      ND1,      ND2,      ND3,      NST,      SWITCH
      INTEGER      IERPLC,   IERROR
C
      REAL*8       FREAL   (M0DIM,N0DIM),FIMAG   (M0DIM,N0DIM)
      REAL*8       PHASER  (N0DIM),      PHASEI  (N0DIM)
C
C
      REAL*8,      DIMENSION(:,:,:),ALLOCATABLE:: WORKINP
      REAL*8,      DIMENSION(:,:,:),ALLOCATABLE:: WORKOUT
      REAL*8,      DIMENSION(:,:),  ALLOCATABLE:: CST
      REAL*8,      DIMENSION(:,:),  ALLOCATABLE:: SNT
C
C
C
C 1.0 INITIALIZATION
C
      SUBNAM    = 'FFTCAL6F'
      IERSUB    = ''
      IERPLC    = 0
      IERROR    = 0
C
      ND3P      = (ND3-1)/2
      ND3Q      =  ND3/2
      ND3QP     =  ND3Q  + 1
      N3PARITY  =  ND3   - 2*(ND3 /2)
C
      IF(ND2 .NE. 6) THEN
         IERSUB   =  SUBNAM
         IERPLC   = -1
         IERROR   =  1
         RETURN
      ENDIF
C
      IF(IABS(SWITCH) .NE. +1) THEN
         IERSUB   =  SUBNAM
         IERPLC   = -2
         IERROR   =  SWITCH
         RETURN
      ENDIF
C
      ROOT3D2   = 0.866025403784438646763723170752936
C
C
C
C 2.0 STORE THE INPUT DATA IN DIFFERENT SIZE WORK ARRAYS
C
C 2.1 CHECK FOR CONSISTENT DIMENSIONS
C
      NDTOTL0   = M0*N0
      NDTOTL1   = ND1*ND2*ND3
      IF(NDTOTL0 .NE. NDTOTL1) THEN
         IERSUB   =  SUBNAM
         IERPLC   = +1
         IERROR   = NDTOTL0 - NDTOTL1
         RETURN
      ENDIF
C
      NDTL0     = ND2*ND3
      NDTL1     = N0 - NST + 1
      IF(NDTL0 .NE. NDTL1) THEN
         IERSUB   =  SUBNAM
         IERPLC   = +2
         IERROR   = NDTL0 - NDTL1
         RETURN
      ENDIF
C
      IF(NST .GT. N0) THEN
         IERSUB   =  SUBNAM
         IERPLC   = +3
         IERROR   =  NST
         RETURN
      ENDIF
C
C
C 2.2 ALLOCATE INTERNAL ARRAYS
C
      ALLOCATE (WORKINP(ND1,ND2,ND3))
      ALLOCATE (WORKOUT(ND1,ND3,ND2))
      ALLOCATE (CST    (ND3,ND2))
      ALLOCATE (SNT    (ND3,ND2))
C
C
C 2.3 FILL IN THE INTERNAL ARRAYS
C
C 2.3.1 FILL IN THE WORK ARRAYS WITH DIMENSIONS (ND1,ND2,ND3) AND (ND1,ND3,ND2)
C
      IVP      = 0
      KVI      = 1
      LVI      = 1
      KVO      = 1
      LVO      = 1
      DO J     = 1,N0
      JV       = J
      DO I     = 1,M0
      IV       = I
      IVP      = IVP + 1
      IF(IVP .GT. ND1) THEN
         IVP      = 1
         KVI      = KVI + 1
         LVO      = LVO + 1
         IF(KVI .GT. ND2) THEN
            KVI       = 1
            LVI       = LVI + 1
            IF(LVI .GT. ND3) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +4
               IERROR   = LVI
               RETURN
            ENDIF
         ENDIF
         IF(LVO .GT. ND3) THEN
            LVO       = 1
            KVO       = KVO + 1
            IF(KVO .GT. ND2) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +5
               IERROR   = KVO
               RETURN
            ENDIF
         ENDIF
      ENDIF
C
      IF    (SWITCH .EQ. +1) THEN
         FR1    = FREAL(IV,JV)
      ELSEIF(SWITCH .EQ. -1) THEN
         FR1    = FIMAG(IV,JV)
      ENDIF
      WORKINP(IVP,KVI,LVI)  = FR1
      WORKOUT(IVP,LVO,KVO)  = 0.0
      ENDDO
      ENDDO
C
C 2.3.2 FILL IN THE PHASE ARRAYS FROM K = NST TO N0
C
      LVP      = 0
      KVP      = 1
      DO J     = NST,N0
      JV       = J
      LVP      = LVP + 1
      IF(LVP .GT. ND3) THEN
         LVP      = 1
         KVP      = KVP + 1
         IF(KVP .GT. ND2) THEN
            IERSUB   =  SUBNAM
            IERPLC   = +6
            IERROR   = KVP
            RETURN
         ENDIF
      ENDIF
      CST(LVP,KVP)  = PHASER(JV)
      SNT(LVP,KVP)  = PHASEI(JV)
      ENDDO
C
C
C
C 3.0 KV = 1: SET FIRST ELEMENT OF MIDDLE INDEX
C
      DO IV              = 1,ND1
      ZVAR0              =           WORKINP(IV,3,  1) +
     &                               WORKINP(IV,5,  1)
      ZVAR1R             =           WORKINP(IV,1,  1) - 0.5*ZVAR0
      ZVAR1I             = +ROOT3D2*(WORKINP(IV,3,  1) -
     &                               WORKINP(IV,5,  1))
      ZVAR1              =           WORKINP(IV,1,  1) +     ZVAR0
      ZVAR2              =           WORKINP(IV,6,  1) +
     &                               WORKINP(IV,2,  1)
      ZVAR3R             =           WORKINP(IV,4,  1) - 0.5*ZVAR2
      ZVAR3I             = +ROOT3D2*(WORKINP(IV,6,  1) -
     &                               WORKINP(IV,2,  1))
      ZVAR3              =           WORKINP(IV,4,  1) +     ZVAR2
C
      WORKOUT(IV,1, 1)   = +ZVAR1  + ZVAR3
      WORKOUT(IV,1, 2)   = +ZVAR1R - ZVAR3R
      WORKOUT(IV,1, 3)   = +ZVAR1R + ZVAR3R
      WORKOUT(IV,1, 4)   = +ZVAR1  - ZVAR3
      WORKOUT(IV,1, 5)   = +ZVAR1I + ZVAR3I
      WORKOUT(IV,1, 6)   = -ZVAR1I + ZVAR3I
      ENDDO
C
C
Ca
C 4.0 FOR KV = 2, ND3P + 1
C
C 4.1 FOR ND1 > (ND3-1)/2
C
      IF    (ND1 .GT. ND3P) THEN
         DO IV              = 1,ND1
            DO K               = 1,ND3P
            KV                 = K + 1
            KPV                = ND3 - K + 1
            CSK1               = 1.0
            CSK2               = CST(KV,1)
            CSK3               = CST(KV,2)
            CSK4               = CST(KV,3)
            CSK5               = CST(KV,4)
            CSK6               = CST(KV,5)
            SNK1               = 0.0
            SNK2               = SNT(KV,1)
            SNK3               = SNT(KV,2)
            SNK4               = SNT(KV,3)
            SNK5               = SNT(KV,4)
            SNK6               = SNT(KV,5)
C
            ZVAR0R             = CSK1*WORKINP(IV,1, KV)  -
     &                                          SNK1*WORKINP(IV,1, KPV)
            ZVAR0I             = CSK1*WORKINP(IV,1, KPV) +
     &                                          SNK1*WORKINP(IV,1, KV)
            ZVAR1R             = CSK2*WORKINP(IV,2, KV)  -
     &                                          SNK2*WORKINP(IV,2, KPV)
            ZVAR1I             = CSK2*WORKINP(IV,2, KPV) +
     &                                          SNK2*WORKINP(IV,2, KV)
            ZVAR2R             = CSK3*WORKINP(IV,3, KV)  -
     &                                          SNK3*WORKINP(IV,3, KPV)
            ZVAR2I             = CSK3*WORKINP(IV,3, KPV) +
     &                                          SNK3*WORKINP(IV,3, KV)
            ZVAR3R             = CSK4*WORKINP(IV,4, KV)  -
     &                                          SNK4*WORKINP(IV,4, KPV)
            ZVAR3I             = CSK4*WORKINP(IV,4, KPV) +
     &                                          SNK4*WORKINP(IV,4, KV)
            ZVAR4R             = CSK5*WORKINP(IV,5, KV)  -
     &                                          SNK5*WORKINP(IV,5, KPV)
            ZVAR4I             = CSK5*WORKINP(IV,5, KPV) +
     &                                          SNK5*WORKINP(IV,5, KV)
            ZVAR5R             = CSK6*WORKINP(IV,6, KV)  -
     &                                          SNK6*WORKINP(IV,6, KPV)
            ZVAR5I             = CSK6*WORKINP(IV,6, KPV) +
     &                                          SNK6*WORKINP(IV,6, KV)
C
            ZVAR0R1            =           ZVAR2R   +     ZVAR4R
            ZVAR0I1            =           ZVAR2I   +     ZVAR4I
            ZVAR0R2            =           ZVAR0R   - 0.5*ZVAR0R1
            ZVAR0I2            =           ZVAR0I   - 0.5*ZVAR0I1
            ZVAR0R3            =  ROOT3D2*(ZVAR2R   -     ZVAR4R)
            ZVAR0I3            =  ROOT3D2*(ZVAR2I   -     ZVAR4I)
C
            ZVAR1R1            =           ZVAR0R   +     ZVAR0R1
            ZVAR1I1            =           ZVAR0I   +     ZVAR0I1
            ZVAR1R2            =           ZVAR0R2  +     ZVAR0I3
            ZVAR1I2            =           ZVAR0I2  -     ZVAR0R3
            ZVAR1R3            =           ZVAR0R2  -     ZVAR0I3
            ZVAR1I3            =           ZVAR0I2  +     ZVAR0R3
C
            ZVAR2R1            =           ZVAR5R   +     ZVAR1R
            ZVAR2I1            =           ZVAR5I   +     ZVAR1I
            ZVAR2R2            =           ZVAR3R   - 0.5*ZVAR2R1
            ZVAR2I2            =           ZVAR3I   - 0.5*ZVAR2I1
            ZVAR2R3            =  ROOT3D2*(ZVAR5R   -     ZVAR1R)
            ZVAR2I3            =  ROOT3D2*(ZVAR5I   -     ZVAR1I)
C
            ZVAR3R1            =           ZVAR3R   +     ZVAR2R1
            ZVAR3I1            =           ZVAR3I   +     ZVAR2I1
            ZVAR3R2            =           ZVAR2R2  +     ZVAR2I3
            ZVAR3I2            =           ZVAR2I2  -     ZVAR2R3
            ZVAR3R3            =           ZVAR2R2  -     ZVAR2I3
            ZVAR3I3            =           ZVAR2I2  +     ZVAR2R3
C
            WORKOUT(IV,KV, 1)  =          +ZVAR1R1  +     ZVAR3R1
            WORKOUT(IV,KPV,1)  =          +ZVAR1R3  -     ZVAR3R3
            WORKOUT(IV,KV, 2)  =          +ZVAR1R2  -     ZVAR3R2
            WORKOUT(IV,KPV,2)  =          +ZVAR1R2  +     ZVAR3R2
            WORKOUT(IV,KV, 3)  =          +ZVAR1R3  +     ZVAR3R3
            WORKOUT(IV,KPV,3)  =          +ZVAR1R1  -     ZVAR3R1
            WORKOUT(IV,KV, 4)  =          -ZVAR1I1  +     ZVAR3I1
            WORKOUT(IV,KPV,4)  =          +ZVAR1I3  +     ZVAR3I3
            WORKOUT(IV,KV, 5)  =          -ZVAR1I2  -     ZVAR3I2
            WORKOUT(IV,KPV,5)  =          +ZVAR1I2  -     ZVAR3I2
            WORKOUT(IV,KV, 6)  =          -ZVAR1I3  +     ZVAR3I3
            WORKOUT(IV,KPV,6)  =          +ZVAR1I1  +     ZVAR3I1
            ENDDO
         ENDDO
C
C
C 4.2 FOR ND1 > (ND3-1)/2
C
      ELSEIF(ND1 .LE. ND3P) THEN
         DO K               = 1,ND3P
         KV                 = K   + 1
         KPV                = ND3 - K + 1
         CSK1               = 1.0
         CSK2               = CST(KV,1)
         CSK3               = CST(KV,2)
         CSK4               = CST(KV,3)
         CSK5               = CST(KV,4)
         CSK6               = CST(KV,5)
         SNK1               = 0.0
         SNK2               = SNT(KV,1)
         SNK3               = SNT(KV,2)
         SNK4               = SNT(KV,3)
         SNK5               = SNT(KV,4)
         SNK6               = SNT(KV,5)
            DO IV              = 1,ND1
            ZVAR0R             = CSK1*WORKINP(IV,1, KV)  -
     &                                          SNK1*WORKINP(IV,1, KPV)
            ZVAR0I             = CSK1*WORKINP(IV,1, KPV) +
     &                                          SNK1*WORKINP(IV,1, KV)
            ZVAR1R             = CSK2*WORKINP(IV,2, KV)  -
     &                                          SNK2*WORKINP(IV,2, KPV)
            ZVAR1I             = CSK2*WORKINP(IV,2, KPV) +
     &                                          SNK2*WORKINP(IV,2, KV)
            ZVAR2R             = CSK3*WORKINP(IV,3, KV)  -
     &                                          SNK3*WORKINP(IV,3, KPV)
            ZVAR2I             = CSK3*WORKINP(IV,3, KPV) +
     &                                          SNK3*WORKINP(IV,3, KV)
            ZVAR3R             = CSK4*WORKINP(IV,4, KV)  -
     &                                          SNK4*WORKINP(IV,4, KPV)
            ZVAR3I             = CSK4*WORKINP(IV,4, KPV) +
     &                                          SNK4*WORKINP(IV,4, KV)
            ZVAR4R             = CSK5*WORKINP(IV,5, KV)  -
     &                                          SNK5*WORKINP(IV,5, KPV)
            ZVAR4I             = CSK5*WORKINP(IV,5, KPV) +
     &                                          SNK5*WORKINP(IV,5, KV)
            ZVAR5R             = CSK6*WORKINP(IV,6, KV)  -
     &                                          SNK6*WORKINP(IV,6, KPV)
            ZVAR5I             = CSK6*WORKINP(IV,6, KPV) +
     &                                          SNK6*WORKINP(IV,6, KV)
C
            ZVAR0R1            =           ZVAR2R   +     ZVAR4R
            ZVAR0I1            =           ZVAR2I   +     ZVAR4I
            ZVAR0R2            =           ZVAR0R   - 0.5*ZVAR0R1
            ZVAR0I2            =           ZVAR0I   - 0.5*ZVAR0I1
            ZVAR0R3            =  ROOT3D2*(ZVAR2R   -     ZVAR4R)
            ZVAR0I3            =  ROOT3D2*(ZVAR2I   -     ZVAR4I)
C
            ZVAR1R1            =           ZVAR0R   +     ZVAR0R1
            ZVAR1I1            =           ZVAR0I   +     ZVAR0I1
            ZVAR1R2            =           ZVAR0R2  +     ZVAR0I3
            ZVAR1I2            =           ZVAR0I2  -     ZVAR0R3
            ZVAR1R3            =           ZVAR0R2  -     ZVAR0I3
            ZVAR1I3            =           ZVAR0I2  +     ZVAR0R3
C
            ZVAR2R1            =           ZVAR5R   +     ZVAR1R
            ZVAR2I1            =           ZVAR5I   +     ZVAR1I
            ZVAR2R2            =           ZVAR3R   - 0.5*ZVAR2R1
            ZVAR2I2            =           ZVAR3I   - 0.5*ZVAR2I1
            ZVAR2R3            =  ROOT3D2*(ZVAR5R   -     ZVAR1R)
            ZVAR2I3            =  ROOT3D2*(ZVAR5I   -     ZVAR1I)
C
            ZVAR3R1            =           ZVAR3R   +     ZVAR2R1
            ZVAR3I1            =           ZVAR3I   +     ZVAR2I1
            ZVAR3R2            =           ZVAR2R2  +     ZVAR2I3
            ZVAR3I2            =           ZVAR2I2  -     ZVAR2R3
            ZVAR3R3            =           ZVAR2R2  -     ZVAR2I3
            ZVAR3I3            =           ZVAR2I2  +     ZVAR2R3
C
            WORKOUT(IV,KV, 1)  =          +ZVAR1R1  +     ZVAR3R1
            WORKOUT(IV,KPV,1)  =          +ZVAR1R3  -     ZVAR3R3
            WORKOUT(IV,KV, 2)  =          +ZVAR1R2  -     ZVAR3R2
            WORKOUT(IV,KPV,2)  =          +ZVAR1R2  +     ZVAR3R2
            WORKOUT(IV,KV, 3)  =          +ZVAR1R3  +     ZVAR3R3
            WORKOUT(IV,KPV,3)  =          +ZVAR1R1  -     ZVAR3R1
            WORKOUT(IV,KV, 4)  =          -ZVAR1I1  +     ZVAR3I1
            WORKOUT(IV,KPV,4)  =          +ZVAR1I3  +     ZVAR3I3
            WORKOUT(IV,KV, 5)  =          -ZVAR1I2  -     ZVAR3I2
            WORKOUT(IV,KPV,5)  =          +ZVAR1I2  -     ZVAR3I2
            WORKOUT(IV,KV, 6)  =          -ZVAR1I3  +     ZVAR3I3
            WORKOUT(IV,KPV,6)  =          +ZVAR1I1  +     ZVAR3I1
            ENDDO
         ENDDO
      ENDIF
C
C
C 4.3 SPECIAL CASE FOR K = ND3/2 WITH ND3 EVEN
C
      IF(N3PARITY .EQ. 0) THEN
         DO  IV               = 1,ND1
         ZVAR0                =           WORKINP(IV,1, ND3QP)
         ZVAR1                =           WORKINP(IV,4, ND3QP)
         ZVAR2                =           WORKINP(IV,3, ND3QP) -
     &                                    WORKINP(IV,5, ND3QP)
         ZVAR2R               =           WORKINP(IV,1, ND3QP) +
     &                                                0.5*ZVAR2
         ZVAR2I               = +ROOT3D2*(WORKINP(IV,3, ND3QP) +
     &                                    WORKINP(IV,5, ND3QP))
         ZVAR3                =           WORKINP(IV,2, ND3QP) +
     &                                    WORKINP(IV,6, ND3QP)
         ZVAR3R               =          -WORKINP(IV,4, ND3QP) -
     &                                                0.5*ZVAR3
         ZVAR3I               = +ROOT3D2*(WORKINP(IV,2, ND3QP) -
     &                                    WORKINP(IV,6, ND3QP))
C
         WORKOUT(IV,ND3QP, 1) = +ZVAR2R + ZVAR3I
         WORKOUT(IV,ND3QP, 2) = +ZVAR0  - ZVAR2
         WORKOUT(IV,ND3QP, 3) = +ZVAR2R - ZVAR3I
         WORKOUT(IV,ND3QP, 4) = +ZVAR3R + ZVAR2I
         WORKOUT(IV,ND3QP, 5) = +ZVAR1  - ZVAR3
         WORKOUT(IV,ND3QP, 6) = +ZVAR3R - ZVAR2I
         ENDDO
      ENDIF
C
C
C
C 5.0 SAVE THE RESULTS AND CLEAN UP
C
C 5.1 SAVE THE RESULTS
C
      IVP      = 0
      KVI      = 1
      LVI      = 1
      KVO      = 1
      LVO      = 1
      DO J     = 1,N0
      JV       = J
      DO I     = 1,M0
      IV       = I
      IVP      = IVP + 1
      IF(IVP .GT. ND1) THEN
         IVP      = 1
         KVI      = KVI + 1
         LVO      = LVO + 1
         IF(KVI .GT. ND2) THEN
            KVI       = 1
            LVI       = LVI + 1
            IF(LVI .GT. ND3) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +7
               IERROR   = LVI
               RETURN
            ENDIF
         ENDIF
         IF(LVO .GT. ND3) THEN
            LVO       = 1
            KVO       = KVO + 1
            IF(KVO .GT. ND2) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +8
               IERROR   = KVO
               RETURN
            ENDIF
         ENDIF
      ENDIF
C
      FR1    = WORKINP(IVP,KVI,LVI)
      FR2    = WORKOUT(IVP,LVO,KVO)
      IF    (SWITCH .EQ. +1) THEN
         FREAL(IV,JV)  = FR1
         FIMAG(IV,JV)  = FR2
      ELSEIF(SWITCH .EQ. -1) THEN
         FREAL(IV,JV)  = FR2
         FIMAG(IV,JV)  = FR1
      ENDIF
      ENDDO
      ENDDO
C
C
C 5.2 DEALLOCATE INTERNAL ARRAYS
C
      DEALLOCATE(WORKINP)
      DEALLOCATE(WORKOUT)
      DEALLOCATE(CST)
      DEALLOCATE(SNT)
C
C
C
C 6.0 RETURN AND END
C
      RETURN
      END
      SUBROUTINE FFTCAL6I(FREAL,FIMAG,PHASER,PHASEI,M0DIM,M0,N0DIM,N0
     &                                      ,ND1,ND2,ND3,NST,SWITCH
     &                                      ,IERSUB,IERPLC,IERROR)
C
C
C ------------------------------------------------------------------------------------
C
C CALCULATE INVERSE FAST FOURIER TRANSFORM OF FREAL AND STORE IN FIMAG
C COOLEY-TOOKEY ALGORITHM FOR SPECIAL CASE WITH ND2 = 6
C
C ----------------------------------------------------------------------------------------
C
      IMPLICIT NONE
C
      CHARACTER*16 SUBNAM
      INTEGER      I,        IV,       K,        KV,       KPV
      INTEGER      IVP,      KVP,      LVP,      J,        JV
      INTEGER      KVI,      KVO,      LVI,      LVO
      INTEGER      ND3P,     ND3Q,     ND3QP
      INTEGER      NDTL0,    NDTL1,    NDTOTL0,  NDTOTL1
      INTEGER      N3PARITY
      REAL*8       ROOT3D2
      REAL*8       FR1,      FR2
      REAL*8       ZVAR0,    ZVAR1,    ZVAR2,    ZVAR3
      REAL*8       ZVAR0R,   ZVAR0I,   ZVAR1R,   ZVAR1I,   ZVAR2R,
     &             ZVAR2I,   ZVAR3R,   ZVAR3I,   ZVAR4R,   ZVAR4I,
     &             ZVAR5R,   ZVAR5I,   ZVAR6R,   ZVAR6I
      REAL*8       ZVAR0R1,  ZVAR0I1,  ZVAR0R2,  ZVAR0I2,
     &             ZVAR0R3,  ZVAR0I3
      REAL*8       ZVAR1R1,  ZVAR1I1,  ZVAR1R2,  ZVAR1I2,
     &             ZVAR1R3,  ZVAR1I3
      REAL*8       ZVAR2R1,  ZVAR2I1,  ZVAR2R2,  ZVAR2I2,
     &             ZVAR2R3,  ZVAR2I3
      REAL*8       ZVAR3R1,  ZVAR3I1,  ZVAR3R2,  ZVAR3I2,
     &             ZVAR3R3,  ZVAR3I3
C
      REAL*8       CSK1,     CSK2,     CSK3,     CSK4,
     &             CSK5,     CSK6
      REAL*8       SNK1,     SNK2,     SNK3,     SNK4,
     &             SNK5,     SNK6
C
      CHARACTER*16 IERSUB
      INTEGER      M0,       N0
      INTEGER      M0DIM,    N0DIM
      INTEGER      ND1,      ND2,      ND3,      NST,      SWITCH
      INTEGER      IERPLC,   IERROR
C
      REAL*8       FREAL   (M0DIM,N0DIM),FIMAG   (M0DIM,N0DIM)
      REAL*8       PHASER  (N0DIM),      PHASEI  (N0DIM)
C
C
      REAL*8,      DIMENSION(:,:,:),ALLOCATABLE:: WORKINP
      REAL*8,      DIMENSION(:,:,:),ALLOCATABLE:: WORKOUT
      REAL*8,      DIMENSION(:,:),  ALLOCATABLE:: CST
      REAL*8,      DIMENSION(:,:),  ALLOCATABLE:: SNT
C
C
C
C 1.0 INITIALIZATION
C
      SUBNAM    = 'FFTCAL6I'
      IERSUB    = ''
      IERPLC    = 0
      IERROR    = 0
C
      ND3P      = (ND3-1)/2
      ND3Q      =  ND3/2
      ND3QP     =  ND3Q  + 1
      N3PARITY  =  ND3   - 2*(ND3 /2)
C
      IF(ND2 .NE. 6) THEN
         IERSUB   =  SUBNAM
         IERPLC   = -1
         IERROR   =  1
         RETURN
      ENDIF
C
      IF(IABS(SWITCH) .NE. +1) THEN
         IERSUB   =  SUBNAM
         IERPLC   = -2
         IERROR   =  SWITCH
         RETURN
      ENDIF
C
      ROOT3D2   = 0.866025403784438646763723170752936
C
C
C
C 2.0 STORE THE INPUT DATA IN DIFFERENT SIZE WORK ARRAYS
C
C 2.1 CHECK FOR CONSISTENT DIMENSIONS
C
      NDTOTL0   = M0*N0
      NDTOTL1   = ND1*ND2*ND3
      IF(NDTOTL0 .NE. NDTOTL1) THEN
         IERSUB   =  SUBNAM
         IERPLC   = +1
         IERROR   = NDTOTL0 - NDTOTL1
         RETURN
      ENDIF
C
      NDTL0     = ND2*ND3
      NDTL1     = N0 - NST + 1
      IF(NDTL0 .NE. NDTL1) THEN
         IERSUB   =  SUBNAM
         IERPLC   = +2
         IERROR   = NDTL0 - NDTL1
         RETURN
      ENDIF
C
      IF(NST .GT. N0) THEN
         IERSUB   =  SUBNAM
         IERPLC   = +3
         IERROR   =  NST
         RETURN
      ENDIF
C
C
C 2.2 ALLOCATE INTERNAL ARRAYS
C
      ALLOCATE (WORKINP(ND1,ND3,ND2))
      ALLOCATE (WORKOUT(ND1,ND2,ND3))
      ALLOCATE (CST    (ND3,ND2))
      ALLOCATE (SNT    (ND3,ND2))
C
C
C 2.3 FILL IN THE INTERNAL ARRAYS
C
C 2.3.1 FILL IN THE WORK ARRAYS WITH DIMENSIONS (ND1,ND3,ND2) AND (ND1,ND2,ND3)
C
      IVP      = 0
      KVI      = 1
      LVI      = 1
      KVO      = 1
      LVO      = 1
      DO J     = 1,N0
      JV       = J
      DO I     = 1,M0
      IV       = I
      IVP      = IVP + 1
      IF(IVP .GT. ND1) THEN
         IVP      = 1
         LVI      = LVI + 1
         KVO      = KVO + 1
         IF(LVI .GT. ND3) THEN
            LVI       = 1
            KVI       = KVI + 1
            IF(KVI .GT. ND2) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +4
               IERROR   = KVI
               RETURN
            ENDIF
         ENDIF
         IF(KVO .GT. ND2) THEN
            KVO       = 1
            LVO       = LVO + 1
            IF(LVO .GT. ND3) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +5
               IERROR   = LVO
               RETURN
            ENDIF
         ENDIF
      ENDIF
C
      IF    (SWITCH .EQ. +1) THEN
         FR1    = FREAL(IV,JV)
      ELSEIF(SWITCH .EQ. -1) THEN
         FR1    = FIMAG(IV,JV)
      ENDIF
      WORKINP(IVP,LVI,KVI)  = FR1
      WORKOUT(IVP,KVO,LVO)  = 0.0
      ENDDO
      ENDDO
C
C 2.3.2 FILL IN THE PHASE ARRAYS FROM K = NST TO N0
C
      LVP      = 0
      KVP      = 1
      DO J     = NST,N0
      JV       = J
      LVP      = LVP + 1
      IF(LVP .GT. ND3) THEN
         LVP      = 1
         KVP      = KVP + 1
         IF(KVP .GT. ND2) THEN
            IERSUB   =  SUBNAM
            IERPLC   = +6
            IERROR   = KVP
            RETURN
         ENDIF
      ENDIF
      CST(LVP,KVP)  = PHASER(JV)
      SNT(LVP,KVP)  = PHASEI(JV)
      ENDDO
C
C
C
C 3.0 KV = 1: SET FIRST ELEMENT OF MIDDLE INDEX
C
      DO IV              = 1,ND1
      ZVAR0              =           WORKINP(IV,1,  3)
      ZVAR0R             =           WORKINP(IV,1,  1) - 0.5*ZVAR0
      ZVAR0I             =  +ROOT3D2*WORKINP(IV,1,  5)
      ZVAR1              =           WORKINP(IV,1,  1) +     ZVAR0
      ZVAR1R             =           ZVAR0R            +     ZVAR0I
      ZVAR1I             =           ZVAR0R            -     ZVAR0I
      ZVAR2              =           WORKINP(IV,1,  2)
      ZVAR2R             =           WORKINP(IV,1,  4) - 0.5*ZVAR2
      ZVAR2I             =  -ROOT3D2*WORKINP(IV,1,  6)
      ZVAR3              =           WORKINP(IV,1,  4) +     ZVAR2
      ZVAR3R             =           ZVAR2R            +     ZVAR2I
      ZVAR3I             =           ZVAR2R            -     ZVAR2I
C
      WORKOUT(IV,1, 1)   = +ZVAR1   + ZVAR3
      WORKOUT(IV,2, 1)   = +ZVAR1R  - ZVAR3R
      WORKOUT(IV,3, 1)   = +ZVAR1I  + ZVAR3I
      WORKOUT(IV,4, 1)   = +ZVAR1   - ZVAR3
      WORKOUT(IV,5, 1)   = +ZVAR1R  + ZVAR3R
      WORKOUT(IV,6, 1)   = +ZVAR1I  - ZVAR3I
      ENDDO
C
C
C
C 4.0 FOR KV = 2, ND3P + 1
C
C 4.1 FOR ND1 > (ND3-1)/2
C
      IF    (ND1 .GT. ND3P) THEN
         DO IV              = 1,ND1
            DO K               = 1,ND3P
            KV                 = K   + 1
            KPV                = ND3 - K + 1
            CSK1               = 1.0
            CSK2               = CST(KV,1)
            CSK3               = CST(KV,2)
            CSK4               = CST(KV,3)
            CSK5               = CST(KV,4)
            CSK6               = CST(KV,5)
            SNK1               = 0.0
            SNK2               = SNT(KV,1)
            SNK3               = SNT(KV,2)
            SNK4               = SNT(KV,3)
            SNK5               = SNT(KV,4)
            SNK6               = SNT(KV,5)
C
            ZVAR1R             = WORKINP(IV,KV,  1)
            ZVAR1I             = WORKINP(IV,KPV, 1)
            ZVAR2R             = WORKINP(IV,KV,  2)
            ZVAR2I             = WORKINP(IV,KPV, 2)
            ZVAR3R             = WORKINP(IV,KV,  3)
            ZVAR3I             = WORKINP(IV,KPV, 3)
            ZVAR4R             = WORKINP(IV,KV,  4)
            ZVAR4I             = WORKINP(IV,KPV, 4)
            ZVAR5R             = WORKINP(IV,KV,  5)
            ZVAR5I             = WORKINP(IV,KPV, 5)
            ZVAR6R             = WORKINP(IV,KV,  6)
            ZVAR6I             = WORKINP(IV,KPV, 6)
C
            ZVAR0R1            =          +ZVAR3R   +          ZVAR2I
            ZVAR0I1            =          -ZVAR5R   +          ZVAR4I
            ZVAR0R2            =          +ZVAR1R   -      0.5*ZVAR0R1
            ZVAR0I2            =          +ZVAR6I   -      0.5*ZVAR0I1
            ZVAR0R3            = +ROOT3D2*(ZVAR3R   -          ZVAR2I)
            ZVAR0I3            = +ROOT3D2*(ZVAR5R   +          ZVAR4I)
C
            ZVAR1R1            =          +ZVAR1R   +          ZVAR0R1
            ZVAR1I1            =          +ZVAR6I   +          ZVAR0I1
            ZVAR1R2            =          +ZVAR0R2  +          ZVAR0I3
            ZVAR1I2            =          +ZVAR0I2  -          ZVAR0R3
            ZVAR1R3            =          +ZVAR0R2  -          ZVAR0I3
            ZVAR1I3            =          +ZVAR0I2  +          ZVAR0R3
C
            ZVAR2R1            =          +ZVAR2R   +          ZVAR1I
            ZVAR2I1            =          -ZVAR6R   +          ZVAR5I
            ZVAR2R2            =          +ZVAR3I   -      0.5*ZVAR2R1
            ZVAR2I2            =          -ZVAR4R   -      0.5*ZVAR2I1
            ZVAR2R3            = -ROOT3D2*(ZVAR2R   -          ZVAR1I)
            ZVAR2I3            = -ROOT3D2*(ZVAR6R   +          ZVAR5I)
C
            ZVAR3R1            =          +ZVAR3I   +          ZVAR2R1
            ZVAR3I1            =          -ZVAR4R   +          ZVAR2I1
            ZVAR3R2            =          +ZVAR2R2  +          ZVAR2I3
            ZVAR3I2            =          +ZVAR2I2  -          ZVAR2R3
            ZVAR3R3            =          +ZVAR2R2  -          ZVAR2I3
            ZVAR3I3            =          +ZVAR2I2  +          ZVAR2R3
C
            ZVAR1R             =          +ZVAR1R1  +          ZVAR3R1
            ZVAR1I             =          +ZVAR1I1  +          ZVAR3I1
            ZVAR2R             =          +ZVAR1R2  -          ZVAR3R2
            ZVAR2I             =          +ZVAR1I2  -          ZVAR3I2
            ZVAR3R             =          +ZVAR1R3  +          ZVAR3R3
            ZVAR3I             =          +ZVAR1I3  +          ZVAR3I3
            ZVAR4R             =          +ZVAR1R1  -          ZVAR3R1
            ZVAR4I             =          +ZVAR1I1  -          ZVAR3I1
            ZVAR5R             =          +ZVAR1R2  +          ZVAR3R2
            ZVAR5I             =          +ZVAR1I2  +          ZVAR3I2
            ZVAR6R             =          +ZVAR1R3  -          ZVAR3R3
            ZVAR6I             =          +ZVAR1I3  -          ZVAR3I3
C
            WORKOUT(IV,1, KV)  =    +CSK1*ZVAR1R    -     SNK1*ZVAR1I
            WORKOUT(IV,1, KPV) =    +CSK1*ZVAR1I    +     SNK1*ZVAR1R
            WORKOUT(IV,2, KV)  =    +CSK2*ZVAR2R    -     SNK2*ZVAR2I
            WORKOUT(IV,2, KPV) =    +CSK2*ZVAR2I    +     SNK2*ZVAR2R
            WORKOUT(IV,3, KV)  =    +CSK3*ZVAR3R    -     SNK3*ZVAR3I
            WORKOUT(IV,3, KPV) =    +CSK3*ZVAR3I    +     SNK3*ZVAR3R
            WORKOUT(IV,4, KV)  =    +CSK4*ZVAR4R    -     SNK4*ZVAR4I
            WORKOUT(IV,4, KPV) =    +CSK4*ZVAR4I    +     SNK4*ZVAR4R
            WORKOUT(IV,5, KV)  =    +CSK5*ZVAR5R    -     SNK5*ZVAR5I
            WORKOUT(IV,5, KPV) =    +CSK5*ZVAR5I    +     SNK5*ZVAR5R
            WORKOUT(IV,6, KV)  =    +CSK6*ZVAR6R    -     SNK6*ZVAR6I
            WORKOUT(IV,6, KPV) =    +CSK6*ZVAR6I    +     SNK6*ZVAR6R
            ENDDO
         ENDDO
C
C
C 4.2 FOR ND1 > (ND3-1)/2
C
      ELSEIF(ND1 .LE. ND3P) THEN
         DO K               = 1,ND3P
         KV                 = K   + 1
         KPV                = ND3 - K + 1
         CSK1               = 1.0
         CSK2               = CST(KV,1)
         CSK3               = CST(KV,2)
         CSK4               = CST(KV,3)
         CSK5               = CST(KV,4)
         CSK6               = CST(KV,5)
         SNK1               = 0.0
         SNK2               = SNT(KV,1)
         SNK3               = SNT(KV,2)
         SNK4               = SNT(KV,3)
         SNK5               = SNT(KV,4)
         SNK6               = SNT(KV,5)
            DO IV              = 1,ND1
            ZVAR1R             = WORKINP(IV,KV,  1)
            ZVAR1I             = WORKINP(IV,KPV, 1)
            ZVAR2R             = WORKINP(IV,KV,  2)
            ZVAR2I             = WORKINP(IV,KPV, 2)
            ZVAR3R             = WORKINP(IV,KV,  3)
            ZVAR3I             = WORKINP(IV,KPV, 3)
            ZVAR4R             = WORKINP(IV,KV,  4)
            ZVAR4I             = WORKINP(IV,KPV, 4)
            ZVAR5R             = WORKINP(IV,KV,  5)
            ZVAR5I             = WORKINP(IV,KPV, 5)
            ZVAR6R             = WORKINP(IV,KV,  6)
            ZVAR6I             = WORKINP(IV,KPV, 6)
C
            ZVAR0R1            =          +ZVAR3R   +          ZVAR2I
            ZVAR0I1            =          -ZVAR5R   +          ZVAR4I
            ZVAR0R2            =          +ZVAR1R   -      0.5*ZVAR0R1
            ZVAR0I2            =          +ZVAR6I   -      0.5*ZVAR0I1
            ZVAR0R3            = +ROOT3D2*(ZVAR3R   -          ZVAR2I)
            ZVAR0I3            = +ROOT3D2*(ZVAR5R   +          ZVAR4I)
C
            ZVAR1R1            =          +ZVAR1R   +          ZVAR0R1
            ZVAR1I1            =          +ZVAR6I   +          ZVAR0I1
            ZVAR1R2            =          +ZVAR0R2  +          ZVAR0I3
            ZVAR1I2            =          +ZVAR0I2  -          ZVAR0R3
            ZVAR1R3            =          +ZVAR0R2  -          ZVAR0I3
            ZVAR1I3            =          +ZVAR0I2  +          ZVAR0R3
C
            ZVAR2R1            =          +ZVAR2R   +          ZVAR1I
            ZVAR2I1            =          -ZVAR6R   +          ZVAR5I
            ZVAR2R2            =          +ZVAR3I   -      0.5*ZVAR2R1
            ZVAR2I2            =          -ZVAR4R   -      0.5*ZVAR2I1
            ZVAR2R3            = -ROOT3D2*(ZVAR2R   -          ZVAR1I)
            ZVAR2I3            = -ROOT3D2*(ZVAR6R   +          ZVAR5I)
C
            ZVAR3R1            =          +ZVAR3I   +          ZVAR2R1
            ZVAR3I1            =          -ZVAR4R   +          ZVAR2I1
            ZVAR3R2            =          +ZVAR2R2  +          ZVAR2I3
            ZVAR3I2            =          +ZVAR2I2  -          ZVAR2R3
            ZVAR3R3            =          +ZVAR2R2  -          ZVAR2I3
            ZVAR3I3            =          +ZVAR2I2  +          ZVAR2R3
C
            ZVAR1R             =          +ZVAR1R1  +          ZVAR3R1
            ZVAR1I             =          +ZVAR1I1  +          ZVAR3I1
            ZVAR2R             =          +ZVAR1R2  -          ZVAR3R2
            ZVAR2I             =          +ZVAR1I2  -          ZVAR3I2
            ZVAR3R             =          +ZVAR1R3  +          ZVAR3R3
            ZVAR3I             =          +ZVAR1I3  +          ZVAR3I3
            ZVAR4R             =          +ZVAR1R1  -          ZVAR3R1
            ZVAR4I             =          +ZVAR1I1  -          ZVAR3I1
            ZVAR5R             =          +ZVAR1R2  +          ZVAR3R2
            ZVAR5I             =          +ZVAR1I2  +          ZVAR3I2
            ZVAR6R             =          +ZVAR1R3  -          ZVAR3R3
            ZVAR6I             =          +ZVAR1I3  -          ZVAR3I3
C
            WORKOUT(IV,1, KV)  =    +CSK1*ZVAR1R    -     SNK1*ZVAR1I
            WORKOUT(IV,1, KPV) =    +CSK1*ZVAR1I    +     SNK1*ZVAR1R
            WORKOUT(IV,2, KV)  =    +CSK2*ZVAR2R    -     SNK2*ZVAR2I
            WORKOUT(IV,2, KPV) =    +CSK2*ZVAR2I    +     SNK2*ZVAR2R
            WORKOUT(IV,3, KV)  =    +CSK3*ZVAR3R    -     SNK3*ZVAR3I
            WORKOUT(IV,3, KPV) =    +CSK3*ZVAR3I    +     SNK3*ZVAR3R
            WORKOUT(IV,4, KV)  =    +CSK4*ZVAR4R    -     SNK4*ZVAR4I
            WORKOUT(IV,4, KPV) =    +CSK4*ZVAR4I    +     SNK4*ZVAR4R
            WORKOUT(IV,5, KV)  =    +CSK5*ZVAR5R    -     SNK5*ZVAR5I
            WORKOUT(IV,5, KPV) =    +CSK5*ZVAR5I    +     SNK5*ZVAR5R
            WORKOUT(IV,6, KV)  =    +CSK6*ZVAR6R    -     SNK6*ZVAR6I
            WORKOUT(IV,6, KPV) =    +CSK6*ZVAR6I    +     SNK6*ZVAR6R
            ENDDO
         ENDDO
      ENDIF
C
C
C 4.3 SPECIAL CASE FOR K = ND3/2 WITH ND3 EVEN
C
      IF(N3PARITY .EQ. 0) THEN
         DO  IV               = 1,ND1
         ZVAR0                =           WORKINP(IV,ND3QP, 2)
         ZVAR1                =           WORKINP(IV,ND3QP, 5)
         ZVAR2                =           WORKINP(IV,ND3QP, 1)  +
     &                                    WORKINP(IV,ND3QP, 3)
         ZVAR2R               =           WORKINP(IV,ND3QP, 2)  -
     &                                                     0.5*ZVAR2
         ZVAR2I               = +ROOT3D2*(WORKINP(IV,ND3QP, 1)  -
     &                                    WORKINP(IV,ND3QP, 3))
         ZVAR3                =           WORKINP(IV,ND3QP, 4)  +
     &                                    WORKINP(IV,ND3QP, 6)
         ZVAR3R               =           WORKINP(IV,ND3QP, 5)  +
     &                                                     0.5*ZVAR3
         ZVAR3I               = -ROOT3D2*(WORKINP(IV,ND3QP, 4)  -
     &                                    WORKINP(IV,ND3QP, 6))
C
         WORKOUT(IV,1, ND3QP) = +ZVAR0  + ZVAR2
         WORKOUT(IV,2, ND3QP) = +ZVAR3R + ZVAR2I
         WORKOUT(IV,3, ND3QP) = -ZVAR2R + ZVAR3I
         WORKOUT(IV,4, ND3QP) = -ZVAR1  + ZVAR3
         WORKOUT(IV,5, ND3QP) = +ZVAR2R + ZVAR3I
         WORKOUT(IV,6, ND3QP) = +ZVAR3R - ZVAR2I
         ENDDO
      ENDIF
C
C
C
C 5.0 SAVE THE RESULTS AND CLEAN UP
C
C 5.1 SAVE THE RESULTS
C
      IVP      = 0
      KVI      = 1
      LVI      = 1
      KVO      = 1
      LVO      = 1
      DO J     = 1,N0
      JV       = J
      DO I     = 1,M0
      IV       = I
      IVP      = IVP + 1
      IF(IVP .GT. ND1) THEN
         IVP      = 1
         LVI      = LVI + 1
         KVO      = KVO + 1
         IF(LVI .GT. ND3) THEN
            LVI       = 1
            KVI       = KVI + 1
            IF(KVI .GT. ND2) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +7
               IERROR   = KVI
               RETURN
            ENDIF
         ENDIF
         IF(KVO .GT. ND2) THEN
            KVO       = 1
            LVO       = LVO + 1
            IF(LVO .GT. ND3) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +8
               IERROR   = LVO
               RETURN
            ENDIF
         ENDIF
      ENDIF
C
      FR1    = WORKINP(IVP,LVI,KVI)
      FR2    = WORKOUT(IVP,KVO,LVO)
      IF    (SWITCH .EQ. +1) THEN
         FREAL(IV,JV)  = FR1
         FIMAG(IV,JV)  = FR2
      ELSEIF(SWITCH .EQ. -1) THEN
         FREAL(IV,JV)  = FR2
         FIMAG(IV,JV)  = FR1
      ENDIF
      ENDDO
      ENDDO
C
C
C 5.2 DEALLOCATE INTERNAL ARRAYS
C
      DEALLOCATE(WORKINP)
      DEALLOCATE(WORKOUT)
      DEALLOCATE(CST)
      DEALLOCATE(SNT)
C
C
C
C 6.0 RETURN AND END
C
      RETURN
      END

      SUBROUTINE FFTCAL7F(FREAL,FIMAG,PHASER,PHASEI,M0DIM,M0,N0DIM,N0
     &                                      ,ND1,ND2,ND3,NST,SWITCH
     &                                      ,IERSUB,IERPLC,IERROR)
C
C
C ------------------------------------------------------------------------------------
C
C CALCULATE FAST FOURIER TRANSFORM OF FREAL AND STORE IN FIMAG
C COOLEY-TOOKEY ALGORITHM FOR CASE WITH ND2 > 6
C
C ----------------------------------------------------------------------------------------
C
      IMPLICIT NONE
C
      CHARACTER*16 SUBNAM
      INTEGER      I,        J,        K,        L
      INTEGER      IV,       JV,       KV,       LV,       JW
      INTEGER      JP,       LP,       JPV,      KPV,      LPV
      INTEGER      IVP,      KVP,      LVP
      INTEGER      KVI,      KVO,      LVI,      LVO
      INTEGER      ND2P,     ND3P,     ND2M,     ND3Q,     ND3QP
      INTEGER      ND2PP,    ND2PM,    ND2PPP
      INTEGER      IARG0,    IARG
      INTEGER      NDTL0,    NDTL1,    NDTOTL0,  NDTOTL1
      INTEGER      N2PARITY, N3PARITY
      REAL*8       CSTV,     SNTV
      REAL*8       FR1,      FR2
      REAL*8       CREAL,    CIMAG
      REAL*8       ZVAR1,    ZVAR2,    ZVAR3,    ZVAR4
C
      CHARACTER*16 IERSUB
      INTEGER      M0,       N0
      INTEGER      M0DIM,    N0DIM
      INTEGER      ND1,      ND2,      ND3,      NST,      SWITCH
      INTEGER      IERPLC,   IERROR,   ADDERR
C
      REAL*8       FREAL   (M0DIM,N0DIM),FIMAG   (M0DIM,N0DIM)
      REAL*8       PHASER  (N0DIM),      PHASEI  (N0DIM)
C
C
      REAL*8,      DIMENSION(:,:,:),ALLOCATABLE:: WORKINP
      REAL*8,      DIMENSION(:,:,:),ALLOCATABLE:: WORKOUT
      REAL*8,      DIMENSION(:,:),  ALLOCATABLE:: CST
      REAL*8,      DIMENSION(:,:),  ALLOCATABLE:: SNT
C
C
C
C 1.0 INITIALIZATION
C
      SUBNAM    = 'FFTCAL7F'
      IERSUB    = ''
      IERPLC    = 0
      IERROR    = 0
      ADDERR    = 100
C
      ND2M      =  ND2   - 1
      ND2P      = (ND2-1)/2
      ND3P      = (ND3-1)/2
      ND2PM     =  ND2P  - 1
      ND2PP     =  ND2P  + 1
      ND2PPP    =  ND2PP + 1
      ND3Q      =  ND3/2
      ND3QP     =  ND3Q  + 1
C
      N2PARITY  =  ND2P  - 2*(ND2P/2)
      N3PARITY  =  ND3   - 2*(ND3 /2)
C
      IF(ND2 .LE. 6) THEN
         IERSUB   =  SUBNAM
         IERPLC   = -1
         IERROR   =  1
         RETURN
      ENDIF
C
      IF(IABS(SWITCH) .NE. +1) THEN
         IERSUB   =  SUBNAM
         IERPLC   = -2
         IERROR   =  SWITCH
         RETURN
      ENDIF
C
C
C
C 2.0 STORE THE INPUT DATA IN DIFFERENT SIZE WORK AARAYS
C
C 2.1 CHECK FOR CONSISTENT DIMENSIONS
C
      NDTOTL0   = M0*N0
      NDTOTL1   = ND1*ND2*ND3
      IF(NDTOTL0 .NE. NDTOTL1) THEN
         IERSUB   =  SUBNAM
         IERPLC   = +1
         IERROR   = NDTOTL0 - NDTOTL1
         RETURN
      ENDIF
C
      NDTL0     = ND2*ND3
      NDTL1     = N0 - NST + 1
      IF(NDTL0 .NE. NDTL1) THEN
         IERSUB   =  SUBNAM
         IERPLC   = +2
         IERROR   = NDTL0 - NDTL1
         RETURN
      ENDIF
C
      IF(NST .GT. N0) THEN
         IERSUB   =  SUBNAM
         IERPLC   = +3
         IERROR   =  NST
         RETURN
      ENDIF
C
C
C 2.2 ALLOCATE INTERNAL ARRAYS
C
      ALLOCATE (WORKINP(ND1,ND2,ND3))
      ALLOCATE (WORKOUT(ND1,ND3,ND2))
      ALLOCATE (CST    (ND3,ND2))
      ALLOCATE (SNT    (ND3,ND2))
C
C
C 2.3 FILL IN THE INTERNAL ARRAYS
C
C 2.3.1 FILL IN THE WORK ARRAYS WITH DIMENSIONS (ND1,ND2,ND3) AND (ND1,ND3,ND2)
C
      IVP      = 0
      KVI      = 1
      LVI      = 1
      KVO      = 1
      LVO      = 1
      DO J     = 1,N0
      JV       = J
      DO I     = 1,M0
      IV       = I
      IVP      = IVP + 1
      IF(IVP .GT. ND1) THEN
         IVP      = 1
         KVI      = KVI + 1
         LVO      = LVO + 1
         IF(KVI .GT. ND2) THEN
            KVI       = 1
            LVI       = LVI + 1
            IF(LVI .GT. ND3) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +4
               IERROR   = LVI
               RETURN
            ENDIF
         ENDIF
         IF(LVO .GT. ND3) THEN
            LVO       = 1
            KVO       = KVO + 1
            IF(KVO .GT. ND2) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +5
               IERROR   = KVO
               RETURN
            ENDIF
         ENDIF
      ENDIF
C
      IF    (SWITCH .EQ. +1) THEN
         FR1    = FREAL(IV,JV)
      ELSEIF(SWITCH .EQ. -1) THEN
         FR1    = FIMAG(IV,JV)
      ENDIF
      WORKINP(IVP,KVI,LVI)  = FR1
      WORKOUT(IVP,LVO,KVO)  = 0.0
      ENDDO
      ENDDO
C
C 2.3.2 FILL IN THE PHASE ARRAYS FROM K = NST TO N0
C
      LVP      = 0
      KVP      = 1
      DO J     = NST,N0
      JV       = J
      LVP      = LVP + 1
      IF(LVP .GT. ND3) THEN
         LVP      = 1
         KVP      = KVP + 1
         IF(KVP .GT. ND2) THEN
            IERSUB   =  SUBNAM
            IERPLC   = +6
            IERROR   = KVP
            RETURN
         ENDIF
      ENDIF
      CST(LVP,KVP)  = PHASER(JV)
      SNT(LVP,KVP)  = PHASEI(JV)
      ENDDO
C
C
C
C 3.0 FOR ND1 > ND3/2
C
      IF(ND1 .GT. ND3P) THEN
C
C 3.1 KV = 1
C
C 3.1.1 RESORT THE INPUT ELEMENTS
C
         DO J               = 1,ND2P
         JV                 = J + 1
         JPV                = ND2 - J + 1
         DO  IV             = 1,ND1
         ZVAR1              = +WORKINP(IV,JV, 1)
         ZVAR2              = +WORKINP(IV,JPV,1)
         WORKINP(IV,JV, 1)  = +ZVAR1              + ZVAR2
         WORKINP(IV,JPV,1)  = +ZVAR1              - ZVAR2
         ENDDO
         ENDDO
C
C 3.1.2 SET FIRST ELEMENT OF MIDDLE INDEX
C
C 3.1.2.1 SET ELEMENTS
         DO L    = 1,ND2P
         LV      = L + 1
         LPV     = ND2 - L + 1
            DO IV    = 1,ND1
            WORKOUT(IV,1,LV)   = WORKINP(IV,1,1)
            WORKOUT(IV,1,LPV)  = 0.0
            ENDDO
C
            DO J               = 1,ND2P
            JV                 = J   + 1
            JPV                = ND2 - J + 1
            IARG0              =  L*J
            IARG               =  IARG0 - ND2*(IARG0/ND2)
            IF(IARG .GT. 0) THEN
               CSTV               =  CST(1,IARG)
               SNTV               =  SNT(1,IARG)
            ELSE
               IERSUB             =  SUBNAM
               IERPLC             = -3
               IERROR             =  ADDERR + IARG
               CSTV               =  1.0
               SNTV               =  0.0
            ENDIF
            DO IV              = 1,ND1
            WORKOUT(IV,1,LV)   = WORKOUT(IV,1,  LV)     +
     &                           WORKINP(IV,JV, 1)*CSTV
            WORKOUT(IV,1,LPV)  = WORKOUT(IV,1,  LPV)    +
     &                           WORKINP(IV,JPV,1)*SNTV
            ENDDO
            ENDDO
         ENDDO
C
C 3.1.2.2 SET FIRST ELEMENT OF MIDDLE AND LAST INDEX
         DO IV              = 1,ND1
         WORKOUT(IV,1,1)    = WORKINP(IV,1,1)
         ENDDO
C
         DO J               = 1,ND2P
         JV                 = J + 1
         DO IV              = 1,ND1
         WORKOUT(IV,1,1)    = WORKOUT(IV,1,1) + WORKINP(IV,JV,1)
         ENDDO
         ENDDO
C
C
C 3.2 fOR KV = 2, ND3P + 1
C
C 3.2.1 SET ELEMENTS
C
         DO K                  = 1,ND3P
         KV                    = K   + 1
         KPV                   = ND3 - K + 1
            DO J                  = 1,ND2M
            JV                    = J + 1
            DO IV                 = 1,ND1
            CREAL                 = WORKINP(IV,JV,KV)
            CIMAG                 = WORKINP(IV,JV,KPV)
            WORKINP(IV,JV,KV)     = CREAL*CST(KV,J) - CIMAG*SNT(KV,J)
            WORKINP(IV,JV,KPV)    = CIMAG*CST(KV,J) + CREAL*SNT(KV,J)
            ENDDO
            ENDDO
C
            DO J                  = 1,ND2P
            JV                    = J   + 1
            JPV                   = ND2 - J + 1
            DO IV                 = 1,ND1
            ZVAR1                 = +WORKINP(IV,JV, KV)
            ZVAR2                 = +WORKINP(IV,JV, KPV)
            WORKINP(IV,JV, KV)    = +WORKINP(IV,JPV,KV)  + ZVAR1
            WORKINP(IV,JV, KPV)   = +WORKINP(IV,JPV,KPV) + ZVAR2
            WORKINP(IV,JPV,KV)    = -WORKINP(IV,JPV,KV)  + ZVAR1
            WORKINP(IV,JPV,KPV)   = -WORKINP(IV,JPV,KPV) + ZVAR2
            ENDDO
            ENDDO
C
            DO L                  = 1,ND2P
            LV                    = L   + 1
            LP                    = ND2 - L
            LPV                   = LP  + 1
               DO IV                 = 1,ND1
               WORKOUT(IV,KV,LV)     = +WORKINP(IV,1,KV)
               WORKOUT(IV,KPV,LP)    = +WORKINP(IV,1,KPV)
               WORKOUT(IV,KPV,L)     =  0.0
               WORKOUT(IV,KV,LPV)    =  0.0
               ENDDO
C
               DO J                  = 1,ND2P
               JV                    = J   + 1
               JPV                   = ND2 - J + 1
               IARG0                 =  L*J
               IARG                  =  IARG0 - ND2*(IARG0/ND2)
               IF(IARG .GT. 0) THEN
                  CSTV                 =  CST(1,IARG)
                  SNTV                 =  SNT(1,IARG)
               ELSE
                  IERSUB               =  SUBNAM
                  IERPLC               = -4
                  IERROR               =  ADDERR + IARG
                  CSTV                 =  1.0
                  SNTV                 =  0.0
               ENDIF
               DO IV                 = 1,ND1
               WORKOUT(IV,KV, LV)    =  WORKOUT(IV,KV, LV)       +
     &                                  WORKINP(IV,JV, KV) *CSTV
               WORKOUT(IV,KPV,LP)    =  WORKOUT(IV,KPV,LP)       +
     &                                  WORKINP(IV,JV, KPV)*CSTV
               WORKOUT(IV,KPV,L)     =  WORKOUT(IV,KPV,L)        -
     &                                  WORKINP(IV,JPV,KV) *SNTV
               WORKOUT(IV,KV, LPV)   =  WORKOUT(IV,KV, LPV)      +
     &                                  WORKINP(IV,JPV,KPV)*SNTV
               ENDDO
               ENDDO
            ENDDO
C
C 3.2.2 SET THE END POINTS
C
            DO IV                 = 1,ND1
            WORKOUT(IV,KV, 1)     =  WORKINP(IV,1,KV)
            WORKOUT(IV,KPV,ND2)   =  WORKINP(IV,1,KPV)
            ENDDO
C
            DO J                  = 1,ND2P
            JV                    = J   + 1
            JPV                   = ND2 - J + 1
            DO IV                 = 1,ND1
            WORKOUT(IV,KV, 1)     =  WORKOUT(IV,KV, 1)   +
     &                               WORKINP(IV,JV, KV)
            WORKOUT(IV,KPV,ND2)   =  WORKOUT(IV,KPV,ND2) +
     &                               WORKINP(IV,JV, KPV)
            ENDDO
            ENDDO
C
C 3.2.3 SET THE FINAL ELEMENT VALUES
C
            DO L                  = 1,ND2P
            LV                    = L   + 1
            LP                    = ND2 - L
            LPV                   = LP  + 1
            DO IV                 = 1,ND1
            ZVAR1                 =  WORKOUT(IV,KV, LV) 
            ZVAR2                 =  WORKOUT(IV,KV, LPV)
            ZVAR3                 =  WORKOUT(IV,KPV,L)
            ZVAR4                 =  WORKOUT(IV,KPV,LP)
            WORKOUT(IV,KV, LV)    = +ZVAR1                - ZVAR2
            WORKOUT(IV,KPV,L)     = +ZVAR1                + ZVAR2
            WORKOUT(IV,KPV,LP)    = -ZVAR3                + ZVAR4
            WORKOUT(IV,KV, LPV)   = -ZVAR3                - ZVAR4
            ENDDO
            ENDDO
         ENDDO
C
C
C 3.3 SPECIAL CASE FOR K = ND3QP WHEN ND3 IS EVEN
C
         IF(N3PARITY .EQ. 0) THEN
C
C 3.3.1 SET ELEMENTS
C
            DO J                  = 1,ND2PM,2
            JV                    = J   + 1
            JW                    = JV  + 1
            JP                    = ND2 - J
            JPV                   = ND2 - J + 1
            DO IV                 = 1,ND1
            ZVAR1                 =  WORKINP(IV,JV, ND3QP)
            ZVAR2                 =  WORKINP(IV,JPV,ND3QP)
            WORKINP(IV,JV, ND3QP) = -ZVAR1                 + ZVAR2
            WORKINP(IV,JPV,ND3QP) = -ZVAR1                 - ZVAR2
            ZVAR3                 =  WORKINP(IV,JW, ND3QP)
            ZVAR4                 =  WORKINP(IV,JP, ND3QP)
            WORKINP(IV,JW, ND3QP) = +ZVAR3                 - ZVAR4
            WORKINP(IV,JP, ND3QP) = +ZVAR3                 + ZVAR4
            ENDDO
            ENDDO
C
C 3.3.2 SET SPECIAL MIDPOINTS
C
            IF(N2PARITY .EQ. 1) THEN
               DO IV                    = 1,ND1
               ZVAR1                    =  WORKINP(IV,ND2PP, ND3QP)
               ZVAR2                    =  WORKINP(IV,ND2PPP,ND3QP)
               WORKINP(IV,ND2PP, ND3QP) = -ZVAR1                + ZVAR2
               WORKINP(IV,ND2PPP,ND3QP) = -ZVAR1                - ZVAR2
               ENDDO
            ENDIF
C
C 3.3.3 SET THE END POINTS
C
            DO LV                  = 1,ND2P
            L                      = LV   - 1
            LP                     = ND2  - L
            LPV                    = LP   + 1
               DO IV                  = 1,ND1
               WORKOUT(IV,ND3QP,LV)   =  WORKINP(IV,1,ND3QP)
               WORKOUT(IV,ND3QP,LP)   =  0.0
               ENDDO
C
               DO J                 = 1,ND2P
               JV                   = J   + 1
               JPV                  = ND2 - J + 1
               IARG0                = J*(L+ND2PP)
               IARG                 =  IARG0 - ND2*(IARG0/ND2)
               IF(IARG .GT. 0) THEN
                  CSTV                =  CST(1,IARG)
                  SNTV                =  SNT(1,IARG)
               ELSE
                  IERSUB              =  SUBNAM
                  IERPLC              = -5
                  IERROR              =  ADDERR + IARG
                  CSTV                =  1.0
                  SNTV                =  0.0
               ENDIF
               DO IV                = 1,ND1
               WORKOUT(IV,ND3QP,LV) =  WORKOUT(IV,ND3QP,LV)          +
     &                                 WORKINP(IV,JV,   ND3QP)*CSTV
               WORKOUT(IV,ND3QP,LP) =  WORKOUT(IV,ND3QP,LP)          +
     &                                 WORKINP(IV,JPV,  ND3QP)*SNTV
               ENDDO
               ENDDO
            ENDDO
C
C 3.3.4 SET THE FINAL ENDPOINT ELEMENT VALUES
C
            DO IV                   = 1,ND1
            WORKOUT(IV,ND3QP,ND2PP) =  WORKINP(IV,1,ND3QP)
            ENDDO
C
            DO J                    = 1,ND2P
            JV                      = J  + 1
            DO IV                   = 1,ND1
            WORKOUT(IV,ND3QP,ND2PP) =  WORKOUT(IV,ND3QP,ND2PP)  +
     &                                 WORKINP(IV,JV,   ND3QP)
            ENDDO
            ENDDO
         ENDIF
C
C
C
C 4.0 FOR ND1 < ND3/2
C
      ELSEIF(ND1 .LE. ND3P) THEN
C
         DO IV                = 1,ND1
C
C 4.1 FOR KV = 1
C
            DO J                 = 1,ND2P
            JV                   = J   + 1
            JPV                  = ND2 - J + 1
            ZVAR1                =  WORKINP(IV,JV, 1)
            ZVAR2                =  WORKINP(IV,JPV,1)
            WORKINP(IV,JV,1)     =  ZVAR1              + ZVAR2
            WORKINP(IV,JPV,1)    =  ZVAR1              - ZVAR2
            ENDDO
C
            DO L                 = 1,ND2P
            LV                   = L   + 1
            LPV                  = ND2 - L + 1
            WORKOUT(IV,1,LV)     =  WORKINP(IV,1,1)
            WORKOUT(IV,1,LPV)    =  0.0
               DO J                 = 1,ND2P
               JV                   = J   + 1
               JPV                  = ND2 - J + 1
               IARG0                =  L*J
               IARG                 =  IARG0 - ND2*(IARG0/ND2)
               IF(IARG .GT. 0) THEN
                  CSTV                =  CST(1,IARG)
                  SNTV                =  SNT(1,IARG)
               ELSE
                  IERSUB              =  SUBNAM
                  IERPLC              = -6
                  IERROR              =  ADDERR + IARG
                  CSTV                =  1.0
                  SNTV                =  0.0
               ENDIF
               WORKOUT(IV,1,LV)     =  WORKOUT(IV,1,  LV)      +
     &                                 WORKINP(IV,JV, 1)*CSTV
               WORKOUT(IV,1,LPV)    =  WORKOUT(IV,1,  LPV)     +
     &                                 WORKINP(IV,JPV,1)*SNTV
               ENDDO             
            ENDDO
C
            WORKOUT(IV,1,1)       =  WORKINP(IV,1,1)
            DO J                  = 1,ND2P
            JV                    = J  + 1
            WORKOUT(IV,1,1)       =  WORKOUT(IV,1,1) + WORKINP(IV,JV,1)
            ENDDO
C
C
C 4.2 FOR K = 2,ND3P+1
C
            DO J                  = 1,ND2M
            JV                    = J  + 1
            DO K                  = 1,ND3P
            KV                    = K   + 1
            KPV                   = ND3 - K + 1
            CREAL                 =  WORKINP(IV,JV,KV)
            CIMAG                 =  WORKINP(IV,JV,KPV)
            WORKINP(IV,JV,KV)     =  CREAL*CST(KV,J)  - CIMAG*SNT(KV,J)
            WORKINP(IV,JV,KPV)    =  CIMAG*CST(KV,J)  + CREAL*SNT(KV,J)
            ENDDO
            ENDDO
C
            DO J                  = 1,ND2P
            JV                    = J   + 1
            JP                    = ND2 - J
            JPV                   = JP  + 1
            DO K                  = 1,ND3P
            KV                    = K   + 1
            KPV                   = ND3 - K + 1
            ZVAR1                 =  WORKINP(IV,JV, KV)
            ZVAR2                 =  WORKINP(IV,JV, KPV)
            WORKINP(IV,JV, KV)    = +WORKINP(IV,JPV,KV)   + ZVAR1
            WORKINP(IV,JV, KPV)   = +WORKINP(IV,JPV,KPV)  + ZVAR2
            WORKINP(IV,JPV,KV)    = -WORKINP(IV,JPV,KV)   + ZVAR1
            WORKINP(IV,JPV,KPV)   = -WORKINP(IV,JPV,KPV)  + ZVAR2
            ENDDO
            ENDDO
C
            DO L                  = 1,ND2P
            LV                    = L   + 1
            LP                    = ND2 - L
            LPV                   = LP  + 1
               DO K                  = 1,ND3P
               KV                    = K    + 1
               KPV                   = ND3  - K + 1
               WORKOUT(IV,KV, LV)    =  WORKINP(IV,1, KV)
               WORKOUT(IV,KPV,LP)    =  WORKINP(IV,1, KPV)
               WORKOUT(IV,KPV,L)     =  0.0
               WORKOUT(IV,KV, LPV)   =  0.0
               ENDDO
C
               DO J                  = 1,ND2P
               JV                    = J   + 1
               JPV                   = ND2 - J + 1
               IARG0                 =  L*J
               IARG                  =  IARG0 - ND2*(IARG0/ND2)
               IF(IARG .GT. 0) THEN
                  CSTV                 =  CST(1,IARG)
                  SNTV                 =  SNT(1,IARG)
               ELSE
                  IERSUB               =  SUBNAM
                  IERPLC               = -7
                  IERROR               =  ADDERR + IARG
                  CSTV                 =  1.0
                  SNTV                 =  0.0
               ENDIF
               DO K                  = 1,ND3P
               KV                    = K   + 1
               KPV                   = ND3 - K + 1
               WORKOUT(IV,KV,LV)     =  WORKOUT(IV,KV, LV)        +
     &                                  WORKINP(IV,JV, KV) *CSTV
               WORKOUT(IV,KPV,LP)    =  WORKOUT(IV,KPV,LP)        +
     &                                  WORKINP(IV,JV, KPV)*CSTV
               WORKOUT(IV,KPV,L)     =  WORKOUT(IV,KPV,L)         -
     &                                  WORKINP(IV,JPV,KV) *SNTV
               WORKOUT(IV,KV ,LPV)   =  WORKOUT(IV,KV ,LPV)       +
     &                                  WORKINP(IV,JPV,KPV)*SNTV
               ENDDO
               ENDDO
            ENDDO
C
            DO K                  = 1,ND3P
            KV                    = K   + 1
            KPV                   = ND3 - K + 1
            WORKOUT(IV,KV, 1)     =  WORKINP(IV,1,KV)
            WORKOUT(IV,KPV,ND2)   =  WORKINP(IV,1,KPV)
            ENDDO
C
            DO J                  = 1,ND2P
            JV                    = J   + 1
            JPV                   = ND2 - J + 1
            DO K                  = 1,ND3P
            KV                    = K   + 1
            KPV                   = ND3 - K + 1
            WORKOUT(IV,KV, 1)     =  WORKOUT(IV,KV, 1)   +
     &                               WORKINP(IV,JV, KV)
            WORKOUT(IV,KPV,ND2)   =  WORKOUT(IV,KPV,ND2) +
     &                               WORKINP(IV,JV, KPV)
            ENDDO
            ENDDO
C
            DO L                  = 1,ND2P
            LV                    = L   + 1
            LP                    = ND2 - L
            LPV                   = LP  + 1
            DO K                  = 1,ND3P
            KV                    = K   + 1
            KPV                   = ND3 - K + 1
            ZVAR1                 =  WORKOUT(IV,KV, LV)
            ZVAR2                 =  WORKOUT(IV,KV, LPV)
            ZVAR3                 =  WORKOUT(IV,KPV,L)
            ZVAR4                 =  WORKOUT(IV,KPV,LP)
            WORKOUT(IV,KV, LV)    = +ZVAR1                - ZVAR2
            WORKOUT(IV,KPV,L)     = +ZVAR1                + ZVAR2
            WORKOUT(IV,KPV,LP)    = -ZVAR3                + ZVAR4
            WORKOUT(IV,KV, LPV)   = -ZVAR3                - ZVAR4
            ENDDO
            ENDDO
C
C
C 4.3 SPECIAL CASE FOR K = ND3QP WHEN ND3 IS EVEN
C
C 4.3.1 SET ELEMENTS
C
            IF(N3PARITY .EQ. 0) THEN
               DO J                  = 1,ND2PM,2
               JV                    = J   + 1
               JW                    = JV  + 1
               JP                    = ND2 - J
               JPV                   = ND2 - J + 1
               ZVAR1                 =  WORKINP(IV,JV, ND3QP)
               ZVAR2                 =  WORKINP(IV,JPV,ND3QP)
               WORKINP(IV,JV, ND3QP) = -ZVAR1                  + ZVAR2
               WORKINP(IV,JPV,ND3QP) = -ZVAR1                  - ZVAR2
               ZVAR3                 = +WORKINP(IV,JW, ND3QP)
               ZVAR4                 = +WORKINP(IV,JP, ND3QP)
               WORKINP(IV,JW, ND3QP) = +ZVAR3                  - ZVAR4
               WORKINP(IV,JP, ND3QP) = +ZVAR3                  + ZVAR4
               ENDDO
C
C 4.3.2 SET THE END POINTS
C
               IF(N2PARITY .EQ. 1) THEN
                  ZVAR1                   =  WORKINP(IV,ND2PP, ND3QP)
                  ZVAR2                   =  WORKINP(IV,ND2PPP,ND3QP)
                  WORKINP(IV,ND2PP, ND3QP)= -ZVAR1               + ZVAR2
                  WORKINP(IV,ND2PPP,ND3QP)= -ZVAR1               - ZVAR2
               ENDIF
C
               DO LV                     = 1,ND2P
               L                         = LV   - 1
               LP                        = ND2P - L
               LPV                       = LP   + 1
               WORKOUT(IV,ND3QP,LV)      =  WORKINP(IV,1,ND3QP)
               WORKOUT(IV,ND3QP,LP)      =  0.0
                  DO J                   = 1,ND2P
                  JV                     = J    + 1
                  JPV                    = ND2  - J + 1
                  IARG0                  = J*(L+ND2PP)
                  IARG                   =  IARG0 - ND2*(IARG0/ND2)
                  IF(IARG .GT. 0) THEN
                     CSTV                  =  CST(1,IARG)
                     SNTV                  =  SNT(1,IARG)
                  ELSE
                     IERSUB                =  SUBNAM
                     IERPLC                = -8
                     IERROR                =  ADDERR + IARG
                     CSTV                  =  1.0
                     SNTV                  =  0.0
                  ENDIF
                  WORKOUT(IV,ND3QP,LV)   = WORKOUT(IV,ND3QP,LV)       +
     &                                     WORKINP(IV,JV,   ND3QP)*CSTV
                  WORKOUT(IV,ND3QP,LP)   = WORKOUT(IV,ND3QP,LP)       + 
     &                                     WORKINP(IV,JPV,  ND3QP)*SNTV
                  ENDDO
               ENDDO
C
C 4.3.3 SET THE FINAL ENDPOINT ELEMENT VALUES
C
               WORKOUT(IV,ND3QP,ND2PP)   =  WORKINP(IV,1,ND3QP)
               DO J                      = 1,ND2P
               JV                        = J  + 1
               WORKOUT(IV,ND3QP,ND2PP)   =  WORKOUT(IV,ND3QP,ND2PP) +
     &                                      WORKINP(IV,JV,   ND3QP)
               ENDDO
            ENDIF
         ENDDO
      ENDIF
C
C
C
C 5.0 SAVE THE RESULTS AND CLEAN UP
C
C 5.1 SAVE THE RESULTS
C
      IVP      = 0
      KVI      = 1
      LVI      = 1
      KVO      = 1
      LVO      = 1
      DO J     = 1,N0
      JV       = J
      DO I     = 1,M0
      IV       = I
      IVP      = IVP + 1
      IF(IVP .GT. ND1) THEN
         IVP      = 1
         KVI      = KVI + 1
         LVO      = LVO + 1
         IF(KVI .GT. ND2) THEN
            KVI       = 1
            LVI       = LVI + 1
            IF(LVI .GT. ND3) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +7
               IERROR   = LVI
               RETURN
            ENDIF
         ENDIF
         IF(LVO .GT. ND3) THEN
            LVO       = 1
            KVO       = KVO + 1
            IF(KVO .GT. ND2) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +8
               IERROR   = KVO
               RETURN
            ENDIF
         ENDIF
      ENDIF
C
      FR1    = WORKINP(IVP,KVI,LVI)
      FR2    = WORKOUT(IVP,LVO,KVO)
      IF    (SWITCH .EQ. +1) THEN
         FREAL(IV,JV)  = FR1
         FIMAG(IV,JV)  = FR2
      ELSEIF(SWITCH .EQ. -1) THEN
         FREAL(IV,JV)  = FR2
         FIMAG(IV,JV)  = FR1
      ENDIF
      ENDDO
      ENDDO
C
C
C 5.2 DEALLOCATE INTERNAL ARRAYS
C
      DEALLOCATE(WORKINP)
      DEALLOCATE(WORKOUT)
      DEALLOCATE(CST)
      DEALLOCATE(SNT)
C
C
C
C 6.0 RETURN AND END
C
      RETURN
      END
      SUBROUTINE FFTCAL7I(FREAL,FIMAG,PHASER,PHASEI,M0DIM,M0,N0DIM,N0
     &                                      ,ND1,ND2,ND3,NST,SWITCH
     &                                      ,IERSUB,IERPLC,IERROR)
C
C
C ------------------------------------------------------------------------------------
C
C CALCULATE INVERSE FAST FOURIER TRANSFORM OF FREAL AND STORE IN FIMAG
C COOLEY-TOOKEY ALGORITHM FOR CASE WITH ND2 > 6
C
C ----------------------------------------------------------------------------------------
C
      IMPLICIT NONE
C
      CHARACTER*16 SUBNAM
      INTEGER      I,        J,        K,        L
      INTEGER      IV,       JV,       KV,       LV,       JW
      INTEGER      JP,       LP,       JPV,      KPV,      LPV
      INTEGER      IVP,      KVP,      LVP
      INTEGER      KVI,      KVO,      LVI,      LVO
      INTEGER      ND2P,     ND3P,     ND2M,     ND3Q,     ND3QP
      INTEGER      ND2PP,    ND2PM,    ND2PPP
      INTEGER      IARG0,    IARG,     ADDERR
      INTEGER      NDTL0,    NDTL1,    NDTOTL0,  NDTOTL1
      INTEGER      N2PARITY, N3PARITY
      REAL*8       CSTV,     SNTV
      REAL*8       FR1,      FR2
      REAL*8       CREAL,    CIMAG
      REAL*8       ZVAR1,    ZVAR2,    ZVAR3,    ZVAR4
      REAL*8       ZVAR3A,   ZVAR4A
C
      CHARACTER*16 IERSUB
      INTEGER      M0,       N0
      INTEGER      M0DIM,    N0DIM
      INTEGER      ND1,      ND2,      ND3,      NST,      SWITCH
      INTEGER      IERPLC,   IERROR
C
      REAL*8       FREAL   (M0DIM,N0DIM),FIMAG   (M0DIM,N0DIM)
      REAL*8       PHASER  (N0DIM),      PHASEI  (N0DIM)
C
C
      REAL*8,      DIMENSION(:,:,:),ALLOCATABLE:: WORKINP
      REAL*8,      DIMENSION(:,:,:),ALLOCATABLE:: WORKOUT
      REAL*8,      DIMENSION(:,:),  ALLOCATABLE:: CST
      REAL*8,      DIMENSION(:,:),  ALLOCATABLE:: SNT
C
C
C
C 1.0 INITIALIZATION
C
      SUBNAM    = 'FFTCAL7I'
      IERSUB    = ''
      IERPLC    = 0
      IERROR    = 0
      ADDERR    = 100
C
      ND2M      =  ND2   - 1
      ND2P      = (ND2-1)/2
      ND3P      = (ND3-1)/2
      ND2PM     =  ND2P  - 1
      ND2PP     =  ND2P  + 1
      ND2PPP    =  ND2PP + 1
      ND3Q      =  ND3/2
      ND3QP     =  ND3Q  + 1
C
      N2PARITY  =  ND2P  - 2*(ND2P/2)
      N3PARITY  =  ND3   - 2*(ND3 /2)
C
      IF(ND2 .LE. 6) THEN
         IERSUB   =  SUBNAM
         IERPLC   = -1
         IERROR   =  1
         RETURN
      ENDIF
C
      IF(IABS(SWITCH) .NE. +1) THEN
         IERSUB   =  SUBNAM
         IERPLC   = -2
         IERROR   =  SWITCH
         RETURN
      ENDIF
C
C
C
C 2.0 STORE THE INPUT DATA IN DIFFERENT SIZE WORK AARAYS
C
C 2.1 CHECK FOR CONSISTENT DIMENSIONS
C
      NDTOTL0   = M0*N0
      NDTOTL1   = ND1*ND2*ND3
      IF(NDTOTL0 .NE. NDTOTL1) THEN
         IERSUB   =  SUBNAM
         IERPLC   = +1
         IERROR   = NDTOTL0 - NDTOTL1
         RETURN
      ENDIF
C
      NDTL0     = ND2*ND3
      NDTL1     = N0 - NST + 1
      IF(NDTL0 .NE. NDTL1) THEN
         IERSUB   =  SUBNAM
         IERPLC   = +2
         IERROR   = NDTL0 - NDTL1
         RETURN
      ENDIF
C
      IF(NST .GT. N0) THEN
         IERSUB   =  SUBNAM
         IERPLC   = +3
         IERROR   =  NST
         RETURN
      ENDIF
C
C
C 2.2 ALLOCATE INTERNAL ARRAYS
C
      ALLOCATE (WORKINP(ND1,ND3,ND2))
      ALLOCATE (WORKOUT(ND1,ND2,ND3))
      ALLOCATE (CST    (ND3,ND2))
      ALLOCATE (SNT    (ND3,ND2))

C
C
C 2.3 FILL IN THE INTERNAL ARRAYS
C
C 2.3.1 FILL IN THE WORK ARRAYS WITH DIMENSIONS (ND1,ND3,ND2) AND (ND1,ND2,ND3)
C
      IVP      = 0
      KVI      = 1
      LVI      = 1
      KVO      = 1
      LVO      = 1
      DO J     = 1,N0
      JV       = J
      DO I     = 1,M0
      IV       = I
      IVP      = IVP + 1
      IF(IVP .GT. ND1) THEN
         IVP      = 1
         LVI      = LVI + 1
         KVO      = KVO + 1
         IF(LVI .GT. ND3) THEN
            LVI       = 1
            KVI       = KVI + 1
            IF(KVI .GT. ND2) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +4
               IERROR   = KVI
               RETURN
            ENDIF
         ENDIF
         IF(KVO .GT. ND2) THEN
            KVO       = 1
            LVO       = LVO + 1
            IF(LVO .GT. ND3) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +5
               IERROR   = LVO
               RETURN
            ENDIF
         ENDIF
      ENDIF
C
      IF    (SWITCH .EQ. +1) THEN
         FR1    = FREAL(IV,JV)
      ELSEIF(SWITCH .EQ. -1) THEN
         FR1    = FIMAG(IV,JV)
      ENDIF
      WORKINP(IVP,LVI,KVI)  = FR1
      WORKOUT(IVP,KVO,LVO)  = 0.0
      ENDDO
      ENDDO
C
C 2.3.2 FILL IN THE PHASE ARRAYS FROM K = NST TO N0
C
      LVP      = 0
      KVP      = 1
      DO J     = NST,N0
      JV       = J
      LVP      = LVP + 1
      IF(LVP .GT. ND3) THEN
         LVP      = 1
         KVP      = KVP + 1
         IF(KVP .GT. ND2) THEN
            IERSUB   =  SUBNAM
            IERPLC   = +6
            IERROR   = KVP
            RETURN
         ENDIF
      ENDIF
      CST(LVP,KVP)  = PHASER(JV)
      SNT(LVP,KVP)  = PHASEI(JV)
      ENDDO
C
C
C
C 3.0 FOR ND1 > ND3/2
C
      IF(ND1 .GT. ND3P) THEN
C
C 3.1 KV = 1
C
C 3.1.1 SET FIRST ELEMENT OF MIDDLE INDEX
C
C 3.1.1.1 SET ELEMENTS
         DO L               = 1,ND2P
         LV                 = L   + 1
         LPV                = ND2 - L + 1
            DO  IV             = 1,ND1
            WORKOUT(IV,LV, 1)  = WORKINP(IV,1,1)
            WORKOUT(IV,LPV,1)  = 0.0
            ENDDO
C
            DO J               = 1,ND2P
            JV                 = J   + 1
            JPV                = ND2 - J + 1
            IARG0              =  L*J
            IARG               =  IARG0 - ND2*(IARG0/ND2)
            IF(IARG .GT. 0) THEN
               CSTV               =  CST(1,IARG)
               SNTV               =  SNT(1,IARG)
            ELSE
               IERSUB             =  SUBNAM
               IERPLC             = -3
               IERROR             =  ADDERR + IARG
               CSTV               =  1.0
               SNTV               =  0.0
            ENDIF
            DO IV              = 1,ND1
            WORKOUT(IV,LV, 1)  = WORKOUT(IV,LV, 1)      +
     &                           WORKINP(IV,1,JV) *CSTV
            WORKOUT(IV,LPV,1)  = WORKOUT(IV,LPV,1)      -
     &                           WORKINP(IV,1,JPV)*SNTV
            ENDDO
            ENDDO
         ENDDO
C
C 3.1.1.2 SET FIRST ELEMENT OF MIDDLE AND LAST INDEX
         DO IV              = 1,ND1
         WORKOUT(IV,1,1)    = WORKINP(IV,1,1)
         ENDDO
C
         DO J               = 1,ND2P
         JV                 = J + 1
         DO IV              = 1,ND1
         WORKOUT(IV,1,1)    = WORKOUT(IV,1,1) + WORKINP(IV,1,JV)
         ENDDO
         ENDDO
C
C 3.1.2 RESORT THE ELEMENTS
C
         DO J               = 1,ND2P
         JV                 = J   + 1
         JPV                = ND2 - J + 1
         DO  IV             = 1,ND1
         ZVAR1              = +WORKOUT(IV,JV, 1)
         ZVAR2              = +WORKOUT(IV,JPV,1)
         WORKOUT(IV,JV, 1)  = +ZVAR1              + ZVAR2
         WORKOUT(IV,JPV,1)  = +ZVAR1              - ZVAR2
         ENDDO
         ENDDO
C
C
C 3.2 FOR KV = 2, ND3P + 1
C
C 3.2.1 SET ELEMENTS
C
         DO K     = 1,ND3P
         KV       = K   + 1
         KPV      = ND3 - K + 1
            DO J                  = 1,ND2P
            JV                    = J   + 1
            JP                    = ND2 - J
            JPV                   = ND2 - J + 1
            DO IV                 = 1,ND1
            ZVAR1                 = +WORKINP(IV,KV, JV)
            ZVAR2                 = +WORKINP(IV,KPV,JP)
            WORKINP(IV,KV, JV)    = +WORKINP(IV,KPV,J)   + ZVAR1
            WORKINP(IV,KPV,JP)    = -WORKINP(IV,KV, JPV) + ZVAR2
            WORKINP(IV,KPV,J)     = -WORKINP(IV,KPV,J)   + ZVAR1
            WORKINP(IV,KV, JPV)   = -WORKINP(IV,KV, JPV) - ZVAR2
            ENDDO
            ENDDO
C
            DO L                  = 1,ND2P
            LV                    = L   + 1
            LP                    = ND2 - L
            LPV                   = ND2 - L + 1
               DO IV                 = 1,ND1
               WORKOUT(IV,LV, KV)    = +WORKINP(IV,KV, 1)
               WORKOUT(IV,LV, KPV)   = +WORKINP(IV,KPV,ND2)
               WORKOUT(IV,LPV,KV)    =  0.0
               WORKOUT(IV,LPV,KPV)   =  0.0
               ENDDO
C
               DO J                  = 1,ND2P
               JV                    = J   + 1
               JP                    = ND2 - J
               JPV                   = ND2 - J + 1
               IARG0                 =  L*J
               IARG                  =  IARG0 - ND2*(IARG0/ND2)
               IF(IARG .GT. 0) THEN
                  CSTV                 =  CST(1,IARG)
                  SNTV                 =  SNT(1,IARG)
               ELSE
                  IERSUB               =  SUBNAM
                  IERPLC               = -4
                  IERROR               =  ADDERR + IARG
                  CSTV                 =  1.0
                  SNTV                 =  0.0
               ENDIF
               DO IV                 = 1,ND1
               WORKOUT(IV,LV, KV)    =  WORKOUT(IV,LV, KV)       +
     &                                  WORKINP(IV,KV, JV) *CSTV
               WORKOUT(IV,LV, KPV)   =  WORKOUT(IV,LV, KPV)      +
     &                                  WORKINP(IV,KPV,JP) *CSTV
               WORKOUT(IV,LPV,KV)    =  WORKOUT(IV,LPV,KV)       +
     &                                  WORKINP(IV,KPV,J)  *SNTV
               WORKOUT(IV,LPV,KPV)   =  WORKOUT(IV,LPV,KPV)      -
     &                                  WORKINP(IV,KV, JPV)*SNTV
               ENDDO
               ENDDO
            ENDDO
C
C 3.2.2 SET THE END POINTS
C
            DO IV                 = 1,ND1
            WORKOUT(IV,1, KV)     =  WORKINP(IV,KV, 1)
            WORKOUT(IV,1, KPV)    =  WORKINP(IV,KPV,ND2)
            ENDDO
C
            DO J                  = 1,ND2P
            JV                    = J   + 1
            JP                    = ND2 - J
            JPV                   = ND2 - J + 1
            DO IV                 = 1,ND1
            WORKOUT(IV,1, KV)     =  WORKOUT(IV,1,  KV)   +
     &                               WORKINP(IV,KV, JV)
            WORKOUT(IV,1, KPV)    =  WORKOUT(IV,1,  KPV)  +
     &                               WORKINP(IV,KPV,JP)
            ENDDO
            ENDDO
C
C 3.2.3 SET THE FINAL ELEMENT VALUES
C
            DO L                  = 1,ND2P
            LV                    = L   + 1
            LP                    = ND2 - L
            LPV                   = ND2 - L + 1
            DO IV                 = 1,ND1
            ZVAR1                 =  WORKOUT(IV,LV, KV) 
            ZVAR2                 =  WORKOUT(IV,LV, KPV)
            ZVAR3                 =  WORKOUT(IV,LPV,KV)
            ZVAR4                 =  WORKOUT(IV,LPV,KPV)
            WORKOUT(IV,LV, KV)    = +ZVAR1              - ZVAR4
            WORKOUT(IV,LV, KPV)   = +ZVAR2              + ZVAR3
            ZVAR3A                =  WORKOUT(IV,LPV,KV) 
            ZVAR4A                =  WORKOUT(IV,LPV,KPV)
            WORKOUT(IV,LPV,KV)    = +ZVAR1              + ZVAR4A
            WORKOUT(IV,LPV,KPV)   = +ZVAR2              - ZVAR3A
            ENDDO
            ENDDO
C
            DO J                  = 1,ND2M
            JV                    = J + 1
            DO IV                 = 1,ND1
            CREAL                 = WORKOUT(IV,JV,KV)
            CIMAG                 = WORKOUT(IV,JV,KPV)
            WORKOUT(IV,JV,KV)     = CREAL*CST(KV,J) - CIMAG*SNT(KV,J)
            WORKOUT(IV,JV,KPV)    = CIMAG*CST(KV,J) + CREAL*SNT(KV,J)
            ENDDO
            ENDDO
         ENDDO
C
C
C 3.3 SPECIAL CASE FOR K = ND3QP WHEN ND3 IS EVEN
C
         IF(N3PARITY .EQ. 0) THEN
C
C 3.3.1 SET THE END POINTS
C
            DO L                  = 1,ND2P
            LV                    = L   + 1
            LPV                   = ND2 - L + 1
            DO IV                 = 1,ND1
            WORKOUT(IV,LV, ND3QP) =  WORKINP(IV,ND3QP,ND2PP)
            WORKOUT(IV,LPV,ND3QP) =  0.0
            ENDDO
            ENDDO
C
            DO L                  = 1,ND2P
            LV                    = L   + 1
            LPV                   = ND2 - L + 1
               DO JV                 = 1,ND2P
               J                     = JV  - 1
               JP                    = ND2 - J
               IARG0                 = L*(J+ND2PP)
               IARG                  =  IARG0 - ND2*(IARG0/ND2)
               IF(IARG .GT. 0) THEN
                  CSTV               =  CST(1,IARG)
                  SNTV               =  SNT(1,IARG)
               ELSE
                  IERSUB             =  SUBNAM
                  IERPLC             = -5
                  IERROR             =  ADDERR + IARG
                  CSTV               =  1.0
                  SNTV               =  0.0
               ENDIF
               DO IV                 = 1,ND1
               WORKOUT(IV,LV, ND3QP) =  WORKOUT(IV,LV,   ND3QP)    +
     &                                  WORKINP(IV,ND3QP,JV) *CSTV
               WORKOUT(IV,LPV,ND3QP) =  WORKOUT(IV,LPV,  ND3QP)    +
     &                                  WORKINP(IV,ND3QP,JP) *SNTV
               ENDDO
               ENDDO
            ENDDO
C
C 3.3.2 SET THE FINAL ENDPOINT ELEMENT VALUES
C
            DO IV                    = 1,ND1
            WORKOUT(IV,1,ND3QP)      =  WORKINP(IV,ND3QP,ND2PP)
            ENDDO
C
            DO JV                    = 1,ND2P
            J                        = JV  - 1
            DO IV                    = 1,ND1
            WORKOUT(IV,1,   ND3QP)   =  WORKOUT(IV,1,    ND3QP)  +
     &                                  WORKINP(IV,ND3QP,JV)
            ENDDO
            ENDDO
C
C 3.3.3 SET ELEMENTS
C
            DO J                   = 1,ND2PM,2
            JV                     = J   + 1
            JW                     = JV  + 1
            JP                     = ND2 - J
            JPV                    = ND2 - J + 1
            DO IV                  = 1,ND1
            ZVAR1                  =  WORKOUT(IV,JV, ND3QP)
            ZVAR2                  =  WORKOUT(IV,JPV,ND3QP)
            WORKOUT(IV,JV, ND3QP)  = -ZVAR1                 + ZVAR2
            WORKOUT(IV,JPV,ND3QP)  = +ZVAR1                 + ZVAR2
            ZVAR3                  =  WORKOUT(IV,JW, ND3QP)
            ZVAR4                  =  WORKOUT(IV,JP, ND3QP)
            WORKOUT(IV,JW, ND3QP)  = +ZVAR3                 - ZVAR4
            WORKOUT(IV,JP, ND3QP)  = -ZVAR3                 - ZVAR4
            ENDDO
            ENDDO
C
C 3.3.4 SET SPECIAL ENDPOINTS
C
            IF(N2PARITY .EQ. 1) THEN
               DO IV                    = 1,ND1
               ZVAR1                    =  WORKOUT(IV,ND2PP, ND3QP)
               ZVAR2                    =  WORKOUT(IV,ND2PPP,ND3QP)
               WORKOUT(IV,ND2PP, ND3QP) = -ZVAR1                 + ZVAR2
               WORKOUT(IV,ND2PPP,ND3QP) = +ZVAR1                 + ZVAR2
               ENDDO
            ENDIF
         ENDIF
C
C
C
C 4.0 FOR ND1 < ND3/2
C
      ELSEIF(ND1 .LE. ND3P) THEN
C
C
         DO IV                = 1,ND1
C
C 4.1 FOR KV = 1
C
C
C 4.1.1 SET FIRST ELEMENT OF MIDDLE INDEX
C
C 4.1.1.1 SET ELEMENTS
            DO L               = 1,ND2P
            LV                 = L   + 1
            LPV                = ND2 - L + 1
            WORKOUT(IV,LV, 1)  = WORKINP(IV,1,1)
            WORKOUT(IV,LPV,1)  = 0.0
C
            DO J               = 1,ND2P
            JV                 = J   + 1
            JPV                = ND2 - J + 1
            IARG0              =  L*J
            IARG               =  IARG0 - ND2*(IARG0/ND2)
            IF(IARG .GT. 0) THEN
               CSTV               =  CST(1,IARG)
               SNTV               =  SNT(1,IARG)
            ELSE
               IERSUB             =  SUBNAM
               IERPLC             = -6
               IERROR             =  ADDERR + IARG
               CSTV               =  1.0
               SNTV               =  0.0
            ENDIF
            WORKOUT(IV,LV, 1)  = WORKOUT(IV,LV, 1)      +
     &                           WORKINP(IV,1,JV) *CSTV
            WORKOUT(IV,LPV,1)  = WORKOUT(IV,LPV,1)      -
     &                           WORKINP(IV,1,JPV)*SNTV
            ENDDO
            ENDDO
C
C 4.1.1.2 SET FIRST ELEMENT OF MIDDLE AND LAST INDEX
            WORKOUT(IV,1,1)    = WORKINP(IV,1,1)
C
            DO J               = 1,ND2P
            JV                 = J + 1
            WORKOUT(IV,1,1)    = WORKOUT(IV,1,1) + WORKINP(IV,1,JV)
            ENDDO
C
C 4.1.2 RESORT THE ELEMENTS
C
            DO J               = 1,ND2P
            JV                 = J   + 1
            JPV                = ND2 - J + 1
            ZVAR1              = +WORKOUT(IV,JV, 1)
            ZVAR2              = +WORKOUT(IV,JPV,1)
            WORKOUT(IV,JV, 1)  = +ZVAR1              + ZVAR2
            WORKOUT(IV,JPV,1)  = +ZVAR1              - ZVAR2
            ENDDO
C
C
C 4.2 FOR KV = 2, ND3P + 1
C
C 4.2.1 SET ELEMENTS
C
            DO J                  = 1,ND2P
            JV                    = J   + 1
            JP                    = ND2 - J
            JPV                   = ND2 - J + 1
            DO K                  = 1,ND3P
            KV                    = K   + 1
            KPV                   = ND3 - K + 1
            ZVAR1                 = +WORKINP(IV,KV, JV)
            ZVAR2                 = +WORKINP(IV,KPV,JP)
            WORKINP(IV,KV, JV)    = +WORKINP(IV,KPV,J)   + ZVAR1
            WORKINP(IV,KPV,JP)    = -WORKINP(IV,KV, JPV) + ZVAR2
            WORKINP(IV,KPV,J)     = -WORKINP(IV,KPV,J)   + ZVAR1
            WORKINP(IV,KV, JPV)   = -WORKINP(IV,KV, JPV) - ZVAR2
            ENDDO
            ENDDO
C
            DO L                  = 1,ND2P
            LV                    = L   + 1
            LP                    = ND2 - L
            LPV                   = LP  + 1
               DO K                  = 1,ND3P
               KV                    = K   + 1
               KPV                   = ND3 - K + 1
               WORKOUT(IV,LV, KV)    = +WORKINP(IV,KV,   1)
               WORKOUT(IV,LV, KPV)   = +WORKINP(IV,KPV,ND2)
               WORKOUT(IV,LPV,KV)    =  0.0
               WORKOUT(IV,LPV,KPV)   =  0.0
               ENDDO
C
               DO J                  = 1,ND2P
               JV                    = J   + 1
               JP                    = ND2 - J
               JPV                   = ND2 - J + 1
               IARG0                 =  L*J
               IARG                  =  IARG0 - ND2*(IARG0/ND2)
               IF(IARG .GT. 0) THEN
                  CSTV               =  CST(1,IARG)
                  SNTV               =  SNT(1,IARG)
               ELSE
                  IERSUB             =  SUBNAM
                  IERPLC             = -7
                  IERROR             =  ADDERR + IARG
                  CSTV               =  1.0
                  SNTV               =  0.0
               ENDIF
               DO K                  = 1,ND3P
               KV                    = K   + 1
               KPV                   = ND3 - K + 1
               WORKOUT(IV,LV, KV)    =  WORKOUT(IV,LV, KV)       +
     &                                  WORKINP(IV,KV, JV) *CSTV
               WORKOUT(IV,LV, KPV)   =  WORKOUT(IV,LV, KPV)      +
     &                                  WORKINP(IV,KPV,JP) *CSTV
               WORKOUT(IV,LPV,KV)    =  WORKOUT(IV,LPV,KV)       +
     &                                  WORKINP(IV,KPV,J)  *SNTV
               WORKOUT(IV,LPV,KPV)   =  WORKOUT(IV,LPV,KPV)      -
     &                                  WORKINP(IV,KV, JPV)*SNTV
               ENDDO
               ENDDO
            ENDDO
C
C 4.2.2 SET THE END POINTS
C
            DO K                  = 1,ND3P
            KV                    = K   + 1
            KPV                   = ND3 - K + 1
            WORKOUT(IV,1, KV)     =  WORKINP(IV,KV, 1)
            WORKOUT(IV,1, KPV)    =  WORKINP(IV,KPV,ND2)
            ENDDO
C
            DO J                  = 1,ND2P
            JV                    = J   + 1
            JP                    = ND2 - J
            JPV                   = ND2 - J + 1
            DO K                  = 1,ND3P
            KV                    = K   + 1
            KPV                   = ND3 - K + 1
            WORKOUT(IV,1, KV)     =  WORKOUT(IV,1,  KV)   +
     &                               WORKINP(IV,KV, JV)
            WORKOUT(IV,1, KPV)    =  WORKOUT(IV,1,  KPV)  +
     &                               WORKINP(IV,KPV,JP)
            ENDDO
            ENDDO
C
C 4.2.3 SET THE FINAL ELEMENT VALUES
C
            DO L                  = 1,ND2P
            LV                    = L   + 1
            LP                    = ND2 - L
            LPV                   = LP  + 1
            DO K                  = 1,ND3P
            KV                    = K   + 1
            KPV                   = ND3 - K + 1
            ZVAR1                 =  WORKOUT(IV,LV, KV) 
            ZVAR2                 =  WORKOUT(IV,LV, KPV)
            ZVAR3                 =  WORKOUT(IV,LPV,KV)
            ZVAR4                 =  WORKOUT(IV,LPV,KPV)
            WORKOUT(IV,LV, KV)    = +ZVAR1              - ZVAR4
            WORKOUT(IV,LV, KPV)   = +ZVAR2              + ZVAR3
            ZVAR3A                =  WORKOUT(IV,LPV,KV)
            ZVAR4A                =  WORKOUT(IV,LPV,KPV)
            WORKOUT(IV,LPV,KV)    = +ZVAR1              + ZVAR4A
            WORKOUT(IV,LPV,KPV)   = +ZVAR2              - ZVAR3A
            ENDDO
            ENDDO
C
            DO J                  = 1,ND2M
            JV                    = J   + 1
            DO K                  = 1,ND3P
            KV                    = K   + 1
            KPV                   = ND3 - K + 1
            CREAL                 = WORKOUT(IV,JV,KV)
            CIMAG                 = WORKOUT(IV,JV,KPV)
            WORKOUT(IV,JV,KV)     = CREAL*CST(KV,J) - CIMAG*SNT(KV,J)
            WORKOUT(IV,JV,KPV)    = CIMAG*CST(KV,J) + CREAL*SNT(KV,J)
            ENDDO
            ENDDO
C
C
C 4.3 SPECIAL CASE FOR K = ND3QP WHEN ND3 IS EVEN
C
            IF(N3PARITY .EQ. 0) THEN
C
C 4.3.1 SET THE END POINTS
C
               DO L                   = 1,ND2P
               LV                     = L   + 1
               LPV                    = ND2 - L + 1
               WORKOUT(IV,LV, ND3QP)  =  WORKINP(IV,ND3QP,ND2PP)
               WORKOUT(IV,LPV,ND3QP)  =  0.0
               ENDDO
C
               DO L                   = 1,ND2P
               LV                     = L   + 1
               LPV                    = ND2 - L + 1
               DO JV                  = 1,ND2P
               J                      = JV  - 1
               JP                     = ND2 - J
               IARG0                  = L*(J+ND2PP)
               IARG                   =  IARG0 - ND2*(IARG0/ND2)
               IF(IARG .GT. 0) THEN
                  CSTV                  =  CST(1,IARG)
                  SNTV                  =  SNT(1,IARG)
               ELSE
                  IERSUB                =  SUBNAM
                  IERPLC                = -8
                  IERROR                =  ADDERR + IARG
                  CSTV                  =  1.0
                  SNTV                  =  0.0
               ENDIF
               WORKOUT(IV,LV, ND3QP)  =  WORKOUT(IV,LV,   ND3QP)     +
     &                                   WORKINP(IV,ND3QP,JV)  *CSTV
               WORKOUT(IV,LPV,ND3QP)  =  WORKOUT(IV,LPV,  ND3QP)     +
     &                                   WORKINP(IV,ND3QP,JP)  *SNTV
               ENDDO
               ENDDO
C
C 4.3.2 SET THE FINAL ENDPOINT ELEMENT VALUES
C
               WORKOUT(IV,1,ND3QP)    =  WORKINP(IV,ND3QP,ND2PP)
C
               DO JV                  = 1,ND2P
               J                      = JV  - 1
               WORKOUT(IV,1,   ND3QP) =  WORKOUT(IV,1,    ND3QP)  +
     &                                   WORKINP(IV,ND3QP,JV)
               ENDDO
C
C 4.3.3 SET ELEMENTS
C
               DO J                   = 1,ND2PM,2
               JV                     = J   + 1
               JW                     = JV  + 1
               JP                     = ND2 - J
               JPV                    = ND2 - J + 1
               ZVAR1                  =  WORKOUT(IV,JV, ND3QP)
               ZVAR2                  =  WORKOUT(IV,JPV,ND3QP)
               WORKOUT(IV,JV, ND3QP)  = -ZVAR1                 + ZVAR2
               WORKOUT(IV,JPV,ND3QP)  = +ZVAR1                 + ZVAR2
               ZVAR3                  =  WORKOUT(IV,JW, ND3QP)
               ZVAR4                  =  WORKOUT(IV,JP, ND3QP)
               WORKOUT(IV,JW, ND3QP)  = +ZVAR3                 - ZVAR4
               WORKOUT(IV,JP, ND3QP)  = -ZVAR3                 - ZVAR4
               ENDDO
C
C 4.3.4 SET SPECIAL ENDPOINTS
C
               IF(N2PARITY .EQ. 1) THEN
                  ZVAR1                   =  WORKOUT(IV,ND2PP, ND3QP)
                  ZVAR2                   =  WORKOUT(IV,ND2PPP,ND3QP)
                  WORKOUT(IV,ND2PP, ND3QP)= -ZVAR1             + ZVAR2
                  WORKOUT(IV,ND2PPP,ND3QP)= +ZVAR1             + ZVAR2
               ENDIF
            ENDIF
         ENDDO
      ENDIF
C
C
C
C 5.0 SAVE THE RESULTS AND CLEAN UP
C
C 5.1 SAVE THE RESULTS
C
      IVP      = 0
      KVI      = 1
      LVI      = 1
      KVO      = 1
      LVO      = 1
      DO J     = 1,N0
      JV       = J
      DO I     = 1,M0
      IV       = I
      IVP      = IVP + 1
      IF(IVP .GT. ND1) THEN
         IVP      = 1
         LVI      = LVI + 1
         KVO      = KVO + 1
         IF(LVI .GT. ND3) THEN
            LVI       = 1
            KVI       = KVI + 1
            IF(KVI .GT. ND2) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +7
               IERROR   = KVI
               RETURN
            ENDIF
         ENDIF
         IF(KVO .GT. ND2) THEN
            KVO       = 1
            LVO       = LVO + 1
            IF(LVO .GT. ND3) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +8
               IERROR   = LVO
               RETURN
            ENDIF
         ENDIF
      ENDIF
C
      FR1    = WORKINP(IVP,LVI,KVI)
      FR2    = WORKOUT(IVP,KVO,LVO)
      IF    (SWITCH .EQ. +1) THEN
         FREAL(IV,JV)  = FR1
         FIMAG(IV,JV)  = FR2
      ELSEIF(SWITCH .EQ. -1) THEN
         FREAL(IV,JV)  = FR2
         FIMAG(IV,JV)  = FR1
      ENDIF
      ENDDO
      ENDDO
C
C
C 5.2 DEALLOCATE INTERNAL ARRAYS
C
      DEALLOCATE(WORKINP)
      DEALLOCATE(WORKOUT)
      DEALLOCATE(CST)
      DEALLOCATE(SNT)
C
C
C
C 6.0 RETURN AND END
C
      RETURN
      END
      SUBROUTINE  PFACTOR(N0,KFACTORS,NKFACTOR,NKF,IERSUB,IERPLC,IERROR)
C
C     FIND FACTORS OF N0 AND STORE IN KFACTORS
C
C     FACTORS OF 6, 4, 2, 3, AND 2 ARE EXTRACTED EXPLICITLY IN THAT ORDER
C     FOR FACTORS OF 5, AND HIGHER PRIMES SWITCHING BETWEEN KJMP = 2 AND 4
C     CATCHES ALL PRIME FACTORS KF BETWEEN 5 AND N0:
C     IT CATCHES ALL ODD NEIGHBOURING PAIRS SKIPPING THOSE
C     THAT SATISFY 6*K + 9  != (2**P)*(3**Q)  FOR SOME P AND Q
C     P MUST BE ZERO SINCE THESE ARE ALL ODD AND IF 6*K + 9 != 3**Q
C     THEN THIS IS FALSE
C
      IMPLICIT NONE
C
      CHARACTER*16 IERSUB
      INTEGER      N0,      NKFACTOR,  NKF
      INTEGER      NPRODCT, NDIFF
      INTEGER      IERPLC,  IERROR
      INTEGER      KFACTORS(NKF)
C
      LOGICAL      FIRST
      CHARACTER*16 SUBNAM
      INTEGER      K,       KF,        KFP,       KJMP
      INTEGER      NP
      INTEGER      KS,      DONE
C
C
C
C 1.0 INITIALIZATION
C
      SUBNAM    = 'PFACTOR'
      IERSUB    = ''
      IERPLC    = 0
      IERROR    = 0
C
      DONE      = 0
      K         = 0
      KF        = 6
      KFP       = 7
      KJMP      = 4
      NP        = N0
C
      FIRST     = .TRUE.
C
C
C
C 2.0 SPECIAL CASE FOR N0 = 1
C
      IF(N0 .EQ. 1) THEN
         NKFACTOR     = 1
         KFACTORS(1)  = 1
         RETURN
      ENDIF
C
C
C
C 3.0 EXTRACT FACTORS OF N0
C
      DO KS     = 1,N0
      IF(DONE .GE. 0) THEN
C
         IF(DONE .EQ. 0) KF   = 6
         IF(DONE .EQ. 1) KF   = 4
         IF(DONE .EQ. 2) KF   = 2
         IF(DONE .EQ. 3) KF   = 3
         IF(DONE .EQ. 4) KF   = 5
         IF(DONE .EQ. 5) KF   = KFP
         IF(DONE .GT. 5) THEN
            IERSUB   =  SUBNAM
            IERPLC   = +1
            IERROR   =  DONE
            NKFACTOR = 0
            RETURN
         ENDIF            
C
C
C 3.2 FACTORS OF 6, 4, 2, 3, 5
C
C 3.2.1 EXTRACT THE FACTOR IF KF IS A FACTOR
C
         IF     (KF .LT. 7) THEN
            IF((KF*(NP/KF) - NP) .EQ. 0) THEN
               K            = K + 1
               KFACTORS(K)  = KF
               NP           = NP/KF
               IF(NP .EQ. 1) DONE  = -1
C
C 3.2.2 INCREMENT TO THE NEXT FACTOR IF KF IS NOT A FACTOR
C
            ELSE
               DONE    = DONE + 1
            ENDIF
C
C
C 3.3 FACTORS OF 5, AND HIGHER PRIMES
C     
         ELSEIF(KF .GE. 7) THEN
C
C 3.3.1 SPECIAL SETUP FOR FIRST TIME THROUGH
C
            IF(FIRST) THEN
               KFP     = 7
               FIRST   = .FALSE.
            ELSE
               KFP     = KF
            ENDIF
C
C 3.3.2 EXTRACT THE FACTOR IF KF IS A FACTOR
C
            IF((KFP*(NP/KFP) - NP) .EQ. 0) THEN
               K            = K + 1
               KFACTORS(K)  = KFP
               NP           = NP/KFP
               IF(NP .EQ. 1) DONE    = -1
C
C 3.3.3 INCREMENT TO THE NEXT PRIME IF KF IS NOT A FACTOR
C
            ELSE
               KFP    = KFP  + KJMP
               KJMP   = 6    - KJMP
            ENDIF
C
C 3.3.4 CHECK LOGIC BEFORE RETURNING TO BEGINNING OF LOOP
C
            IF(KFP .LT. 7) THEN
               IERSUB   =  SUBNAM
               IERPLC   = +2
               IF(KF .EQ. 0) IERROR  = -1
               IF(KF .NE. 0) IERROR  =  KF
               NKFACTOR = 0
               RETURN
            ENDIF
         ENDIF
C
C
C 3.4 EXIT THE LOOP IF ALL FACTORS ARE FOUND
C
      ELSEIF(DONE .LT. 0) THEN
         GO TO 100
      ENDIF
      ENDDO
C
C
C
C 4.0 UPDATE AND CHECK THE FINAL NFACTOR
C
 100  CONTINUE
C
C 4.1 UPDATE THE NUMBER OF FACTORS
C
      NKFACTOR   = K
      IF(NKFACTOR .LE. 0) THEN
         IERSUB  =  SUBNAM
         IERPLC  = +3
         IF(NKFACTOR .EQ. 0) IERROR  =  -1
         IF(NKFACTOR .NE. 0) IERROR  =   NKFACTOR
      ENDIF
C
C
C 4.2 CHECK THE FACTORS
C
      NPRODCT  = 1
      DO K     = 1,NKFACTOR
      NPRODCT  = NPRODCT*KFACTORS(K)
      ENDDO
C
      NDIFF    = NPRODCT - N0
      IF(NDIFF .NE. 0) THEN
         IERSUB  =  SUBNAM
         IERPLC  = +4
         IERROR  = NDIFF
         RETURN
      ENDIF      
C
C
C
C 5.0 RETURN AND END
C
      RETURN
      END
      SUBROUTINE PHASE(PHASER,PHASEI,KFACTR,NKFACT,NKF,N0,NDIM
     &                       ,INITL,KUNIT,IESUBR,INFLAG,IEFLAG)
C
      IMPLICIT NONE
C
      CHARACTER*16 IESUBR
      CHARACTER*16 SUBNAM,  IERSB0
      INTEGER      NCHEK,   NPHASER, NPHASEI
      INTEGER      NKFACT,  NKF,     KFACTR(NKF)
      INTEGER      N0,      NDIM,    INITL
      INTEGER      KUNIT,   INFLAG,  IEFLAG
      INTEGER      IERPL0,  IEROR0
      REAL*8       PHASER(NDIM),     PHASEI(NDIM)
C
      SAVE         NCHEK
      DATA         NCHEK/0/
C
C
C
C 1.0 iNITIALIZATION
C
C 1.1 INITIALIZE ERROR FLAGS
C
      SUBNAM   = 'PHASE'
      IESUBR   = ''
      INFLAG   = 0
      IEFLAG   = 0
C
C
C 1.2 CHECK THE INPUT START FLAG
C     NCHEK ENSURES THE PHASE ARRAYS ARE ALWAYS COMPUTED THE FIRST TIME
C
      IF    (INITL .EQ. +1  .AND.  NCHEK   .EQ.  0) THEN 
         IESUBR   =  SUBNAM
         INFLAG   = +1
         IEFLAG   = -1
         INITL    =  0
      ENDIF
C
      IF(IABS(INITL) .EQ. +1) THEN
         NPHASER  = NINT(PHASER(N0))
         NPHASEI  = NINT(PHASEI(N0))
         IF(NPHASER .NE. N0  .OR.   NPHASEI .NE. N0)  THEN
            IESUBR   =  SUBNAM
            INFLAG   = +2
            IEFLAG   = -2
            RETURN
         ENDIF
      ENDIF
C
C
C
C 2.0 CALCULATE THE INDICES FOR THE PHASES
C
C 2.1 CALCULATE THE FACTORS
C
      CALL PFACTOR(N0,KFACTR,NKFACT,NKF,IERSB0,IERPL0,IEROR0)
C
C
C 2.2 RETURN IF THERE IS AN ERROR
C
C 2.2.1 ERROR IN PFACTOR
C
      IF(IEROR0 .NE. 0) THEN
         IESUBR   =  IERSB0
         INFLAG   =  IERPL0
         IEFLAG   =  IEROR0
         RETURN
      ENDIF
C
C 2.2.2 NUMBER OF FACTORS FOUND IS ZERO
C
      IF(NKFACT .LE. 0) THEN
         IESUBR   =  SUBNAM
         INFLAG   = +3
         IEFLAG   = +1
         RETURN
      ENDIF
C
C
C
C 3.0 CALCULATE THE PHASES ON THE FIRST PASS
C
      IF(INITL .EQ.  0) THEN
         CALL CPHASE(KFACTR,N0,NDIM,NKFACT,NKF,PHASER,PHASEI)
      ENDIF
C
C
C
C 4.0 RESET NCHEK
C
      NCHEK   = +1
C
C
C
C 5.0 RETURN AND END
C
      RETURN
      END
      SUBROUTINE CPHASE(KFACTORS,N0,NDIM,NKFACTOR,NKF,PHASER,PHASEI)
C
      IMPLICIT NONE
C
      REAL*8,   PARAMETER:: PI = 3.1415926535897932385
C
      INTEGER   NP,       LV,       LPV,    NARG
      INTEGER   K,        KFV
      INTEGER   J1,       J1V,      J2,     J2V
      INTEGER   I,        IV
      REAL*8    TWOPI,    ARG
C
      INTEGER   NKFACTOR, NKF
      INTEGER   KFACTORS(NKF)
      INTEGER   N0,       NDIM
      REAL*8    PHASEV
      REAL*8    PHASER(NDIM), PHASEI(NDIM)
C
C
C
C 1.0 INITIALIZATION
C
      NP     = N0
      LV     =  1
      NARG   =  0
C
      TWOPI  = 2.0*PI
      ARG    = TWOPI/DFLOAT(N0)
C
C
C
C 2.0 RUN THROUGH FACTORS OF N0 AND CONSTRUCT PHASES
C
      DO K       = 1,NKFACTOR
C
C
C 2.1 FACTORS OF 2, 3, 4, 5, AND 6
C
         KFV        = KFACTORS(K)
         NP         = NP/KFV
         LPV        = LV
         DO J1V       = 1,KFV-1
         J1           = J1V
         DO J2V       = 1,NP
         J2           = J2V - 1
         NARG         = LV
         PHASER(NARG) = ARG*DFLOAT(J1*J2)
         LV           = LV  + 1
         ENDDO
         ENDDO
C
C
C 2.2 PRIME FACTORS BEYOND 6
C
         IF(KFV .GE. 7) THEN
            LV           = LPV
            DO J1V       = 1,KFV-1
            J1           = J1V
            NARG         = LV
            PHASER(NARG) = ARG*DFLOAT(J1*NP)
            LV           = LV + NP
            ENDDO
         ENDIF
C
C
C 2.3 RESET THE ARGUMENT
C
         ARG     = ARG*KFV
      ENDDO
C
C
C
C 3.0 SET THE SINE AND COSINE OF THE PHASE
C
      DO I         = 1, N0-1
      IV           = I
      PHASEV       =  PHASER(IV)
      PHASEI(IV)   = -SIN(PHASEV)
      PHASER(IV)   = +COS(PHASEV)
      ENDDO
C
C
C     CHECK FOR CONSISTENCY
C
      PHASER(N0) = N0
      PHASEI(N0) = N0
C
C
C
C 4.0 RETURN AND END
C
      RETURN
      END
