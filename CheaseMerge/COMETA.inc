C*COMDECK COMETA
C
***********************************************************************
*                                                                     *
*     C.07         LINEAR STABLITY CODE QUANTITES                     *
*                                                                     *
*     ASPCTR  : ASPECT RATIO OF PLASMA SHAPE = 1 / ASPCT              *
*     EQL     : STORAGE SPACE FOR EQ'S IN DIRECT SPACE                *
*     RM      : REAL POLOIDAL MODE NUMBER                             *
*     RN      : REAL TOROIDAL MODE NUMBER                             *
*                                                                     *
*                       ***************                               *
*                                                                     *
*  METRIC TENSOR QUANTITIES:  (FUNCTIONS OF S AND CHI)                *
*                                                                     *
*            G11L  =  S (1)- S (1) TENSOR ELEMENT ON INTEGER MESH (CS)*
*                                  L MEANS LOW INDICES                *
*            G22L  = CHI(2)-CHI(2) TENSOR ELEMENT ON INTEGER MESH (CS)*
*            G33L  = PHI(3)-PHI(3) TENSOR ELEMENT ON INTEGER MESH (CS)*
*            G12L  =  S (1)-CHI(2) TENSOR ELEMENT ON INTEGER MESH (CS)*
*                                                                     *
*            ....M = ................AS ABOVE.....ON   HALF  MESH (CSM)
*                                                                     *
*           D..... = D*G...    WHERE D IS THE INVERSE OF THE JACOBIAN *
*                                                                     *
*           J..... = J*G...    WHERE J IS THE JACOBIAN ( 1/D )        *
*                                                                     *
*            JACOBI= JACOBIAN ON INTEGER MESH (CS)                    *
*                                                                     *
*            JACOBM= JACOBIAN ON   HALF  MESH (CSM)                   *
*LIUYQ                                                                     *
*            JACOBINV = INVERSE JACOBIAN ON INTEGER MESH (CS)         *
*                                                                     *
*            FRM = FOURIER COEFFICIENTS FOR R                         *
*            FZM = FOURIER COEFFICIENTS FOR Z                         *
*                                                                     *
*                       ***************                               *
*                                                                     *
*      VARIABLES IN X-VECTOR:                                         *
*                                                                     *
*  J2U = COMPONENT OF CONTRAVARIANT (UP) CURRENT DEN.                 *
*  J3U = COMPONENT OF CONTRAVARIANT (UP) CURRENT DEN.                 *
*                                                                     *
*      VARIABLES IN Y-VECTOR:                                         *
*                                                                     *
*  B2U = COMPONENT OF CONTRAVARIANT (UP) MAGN.FIELD                   *
*  B3U = COMPONENT OF CONTRAVARIANT (UP) MAGN.FIELD                   *
*  PRE = PRESSURE                                                     *
*                                                                     *
*  J2U,J3U     ARE DEFINED ON THE INTEGER (CS) MESH                   *
*  B2U,B3U,PRE ARE DEFINED ON THE HALF   (CSM) MESH                   *
*                                                                     *
*  FOR MORE DETAILS SEE SUBROUTINE GIJLIN AND [1], TABLE 3            *
*                                                                     *
***********************************************************************
C
         COMMON /COMETA/
     R   ASPCTR,  EQI,      EQL,     EQ3,
     C   RM,      RN,       DG11L,   DG22L,    
     C   DG33L,   DG12L,    DG11LM,  DG22LM,   
     C   DG33LM,  DG12LM,   JG11L,   JG22L,   
     C   JG33L,   JG12L,    JG11LM,  JG22LM,  
     C   JG33LM,  JG12LM,   JACOBI,  JACOBM,
     C   JACOBINV,FRM,      FZM,   
     C   B2E,     B3E,      B2U,     B3U,    
     C   DPEDS,   DPEDSM,   J2E,     J3E,     
     C   J2U,     J3U,      PEQ,     PRE,
     C   GCHDZ,   GSDZ,     GBZ,     GBR,
     C   GCHDZM,  GSDZM,    GBZM,    GBRM
         COMMON /FLOW/
     C   IDIY2,   IDIY2M,   IDIY3,   IDIY3M,
     C   IG122,   IG122M,   IG123,   IG123M,
     C   INXX,    INXXM,    INXY,    INXYM,
     C   INYY,    INYYM,    INZZ,    INZZM,
     C   IJ0QX,   IJ0QXM,   IJ0QY,   IJ0QYM,
     C   IGPX2,   IGPX2M,   IGPX3,   IGPX3M,
     C   IGPY2,   IGPY2M,   IGPY3,   IGPY3M,
     C   IDRXX,   IDRXXM,   IRXZ,    IRXZM,
     C   IDRYX,   IDRYXM,   IRYX,    IRYXM,
     C   IDRZX,   IDRZXM,   IRZY,    IRZYM,
     C   VISXZ,   VISXZM,   VISYZ,   VISYZM,
     C   IVS11,   IVS11M,   IVS12,   IVS12M,
     C   IVS21,   IVS21M,   IVS22,   IVS22M,
     C   GSFC,    GSFCM,    GSCC,    GSCCM,
     C   GSFS,    GSFSM,    GSCS,    GSCSM,
     C   GCFC,    GCFCM,    GCFS,    GCFSM
C
         COMMON /FLOW1/
     C   EQRHO,   EQRHOM,   DRHOS,   DRHOSM,
     C   EQROT,   EQROTM,   DROT,    DROTM,
     C   FEQ,     FEQM,     IWSQ1,   IWSQ1M,
     C   IWSQ2,   IWSQ2M,   IWSQ3,   IWSQ3M,
     C   IJ0QZ,   IJ0QZM,   JACOF,   JACOFM,
     C   B2F,     B2FM,     B3F,     B3FM,
     C   JACOS,   JACOSM,   IGF22,   IGF22M,
     C   B3FC,    B3FCM,    B2FC,    B2FCM,
     C   DJCOF,   DJCOFM
C
         DIMENSION
     R   EQI(NPMGS*NTP1,32),      EQL(NPMGS*NTP1,22),  
     R   EQ3(NPMGS*NTP1,17),
     R   RM(MPSMAX,NPSMAX),       RN(NPSMAX)
C
         COMPLEX
     C   DG11L(NPISO,MPSMAX),        DG22L(NPISO,MPSMAX),
     C   DG33L(NPISO,MPSMAX),        DG12L(NPISO,MPSMAX),
     C   DG11LM(NPISO,MPSMAX),       DG22LM(NPISO,MPSMAX),
     C   DG33LM(NPISO,MPSMAX),       DG12LM(NPISO,MPSMAX),
     C   JG11L(NPISO,MPSMAX),        JG22L(NPISO,MPSMAX),
     C   JG33L(NPISO,MPSMAX),        JG12L(NPISO,MPSMAX),
     C   JG11LM(NPISO,MPSMAX),       JG22LM(NPISO,MPSMAX),
     C   JG33LM(NPISO,MPSMAX),       JG12LM(NPISO,MPSMAX),
     C   JACOBI(NPISO,MPSMAX),       JACOBM(NPISO,MPSMAX),
     C   JACOBINV(NPISO,MPSMAX),     
     C   FRM(NPISO,MPSMAX),          FZM(NPISO,MPSMAX),
     C   B2E(NPISO,MPSMAX,NPSMAX),   B3E(NPISO,MPSMAX,NPSMAX),
     C   B2U(NPISO,MPSMAX,NPSMAX),   B3U(NPISO,MPSMAX,NPSMAX),
     C   DPEDS(NPISO,MPSMAX,NPSMAX), DPEDSM(NPISO,MPSMAX,NPSMAX),
     C   J2E(NPISO,MPSMAX,NPSMAX),   J3E(NPISO,MPSMAX,NPSMAX),
     C   J2U(NPISO,MPSMAX,NPSMAX),   J3U(NPISO,MPSMAX,NPSMAX),
     C   PEQ(NPISO,MPSMAX,NPSMAX),   PRE(NPISO,MPSMAX,NPSMAX),
     C   GCHDZ(NPISO,MPSMAX),        GSDZ(NPISO,MPSMAX),
     C   GCHDZM(NPISO,MPSMAX),       GSDZM(NPISO,MPSMAX),
     C   GBZ(NPISO,MPSMAX),          GBR(NPISO,MPSMAX),
     C   GBZM(NPISO,MPSMAX),         GBRM(NPISO,MPSMAX)
         COMPLEX
     C   IDIY2(NPISO,MPSMAX),        IDIY2M(NPISO,MPSMAX),
     C   IDIY3(NPISO,MPSMAX),        IDIY3M(NPISO,MPSMAX),
     C   IG122(NPISO,MPSMAX),        IG122M(NPISO,MPSMAX),
     C   IG123(NPISO,MPSMAX),        IG123M(NPISO,MPSMAX),
     C   INXX(NPISO,MPSMAX),         INXXM(NPISO,MPSMAX),
     C   INXY(NPISO,MPSMAX),         INXYM(NPISO,MPSMAX),
     C   INYY(NPISO,MPSMAX),         INYYM(NPISO,MPSMAX),
     C   INZZ(NPISO,MPSMAX),         INZZM(NPISO,MPSMAX),
     C   IJ0QX(NPISO,MPSMAX),        IJ0QXM(NPISO,MPSMAX),
     C   IJ0QY(NPISO,MPSMAX),        IJ0QYM(NPISO,MPSMAX),
     C   IGPX2(NPISO,MPSMAX),        IGPX2M(NPISO,MPSMAX),
     C   IGPX3(NPISO,MPSMAX),        IGPX3M(NPISO,MPSMAX),
     C   IGPY2(NPISO,MPSMAX),        IGPY2M(NPISO,MPSMAX),
     C   IGPY3(NPISO,MPSMAX),        IGPY3M(NPISO,MPSMAX),
     C   IDRXX(NPISO,MPSMAX),        IDRXXM(NPISO,MPSMAX),
     C   IRXZ(NPISO,MPSMAX),         IRXZM(NPISO,MPSMAX),
     C   IDRYX(NPISO,MPSMAX),        IDRYXM(NPISO,MPSMAX)
         COMPLEX
     C   IRYX(NPISO,MPSMAX),         IRYXM(NPISO,MPSMAX),
     C   IDRZX(NPISO,MPSMAX),        IDRZXM(NPISO,MPSMAX),
     C   IRZY(NPISO,MPSMAX),         IRZYM(NPISO,MPSMAX),
     C   VISXZ(NPISO,MPSMAX),        VISXZM(NPISO,MPSMAX),
     C   VISYZ(NPISO,MPSMAX),        VISYZM(NPISO,MPSMAX),
     C   IVS11(NPISO,MPSMAX),        IVS11M(NPISO,MPSMAX),
     C   IVS12(NPISO,MPSMAX),        IVS12M(NPISO,MPSMAX),
     C   IVS21(NPISO,MPSMAX),        IVS21M(NPISO,MPSMAX),
     C   IVS22(NPISO,MPSMAX),        IVS22M(NPISO,MPSMAX),
     C   GSFC(NPISO,MPSMAX),         GSFCM(NPISO,MPSMAX),
     C   GSCC(NPISO,MPSMAX),         GSCCM(NPISO,MPSMAX),
     C   GSFS(NPISO,MPSMAX),         GSFSM(NPISO,MPSMAX),
     C   GSCS(NPISO,MPSMAX),         GSCSM(NPISO,MPSMAX),
     C   GCFC(NPISO,MPSMAX),         GCFCM(NPISO,MPSMAX),
     C   GCFS(NPISO,MPSMAX),         GCFSM(NPISO,MPSMAX)
         COMPLEX
     C   EQRHO(NPISO,MPSMAX),        EQRHOM(NPISO,MPSMAX),
     C   DRHOS(NPISO,MPSMAX),        DRHOSM(NPISO,MPSMAX),
     C   EQROT(NPISO,MPSMAX),        EQROTM(NPISO,MPSMAX),
     C   DROT(NPISO,MPSMAX),         DROTM(NPISO,MPSMAX),
     C   FEQ(NPISO,MPSMAX),          FEQM(NPISO,MPSMAX),
     C   IWSQ1(NPISO,MPSMAX),        IWSQ1M(NPISO,MPSMAX),
     C   IWSQ2(NPISO,MPSMAX),        IWSQ2M(NPISO,MPSMAX),
     C   IWSQ3(NPISO,MPSMAX),        IWSQ3M(NPISO,MPSMAX),
     C   IJ0QZ(NPISO,MPSMAX),        IJ0QZM(NPISO,MPSMAX),
     C   JACOF(NPISO,MPSMAX),        JACOFM(NPISO,MPSMAX),
     C   B2F(NPISO,MPSMAX),          B2FM(NPISO,MPSMAX),
     C   B3F(NPISO,MPSMAX),          B3FM(NPISO,MPSMAX),
     C   JACOS(NPISO,MPSMAX),        JACOSM(NPISO,MPSMAX),
     C   IGF22(NPISO,MPSMAX),        IGF22M(NPISO,MPSMAX),
     C   B3FC(NPISO,MPSMAX),         B3FCM(NPISO,MPSMAX),
     C   B2FC(NPISO,MPSMAX),         B2FCM(NPISO,MPSMAX),
     C   DJCOF(NPISO,MPSMAX),        DJCOFM(NPISO,MPSMAX)
