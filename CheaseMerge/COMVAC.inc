C*COMDECK COMVAC
C
***********************************************************************
*                                                                     *
*     C.18         VACUUM QUANTITIES FOR LINEAR RESISTIVE STABILITY   *
*                                                                     *
*                       ***************                               *
*                                                                     *
*  METRIC TENSOR QUANTITIES:                                          *
*                                                                     *
*       DEFINITIONS : R  = R0W  + S * R1                              *
*                     Z  = RZ0W + S * Z1                              *
*                     R1 = RSURF(CHI) - R0W                           *
*                     Z1 = ZSURF(CHI) - RZ0W                          *
*                                                                     *
*                     WHERE (RSURF(CHI); ZSURF(CHI)) IS THE PLASMA    *
*                     SURFACE                                         *
*                                                                     *
*                     VECT = (R; Z; R * PHI)                          *
*                                                                     *
*                                                                     *
*       QUANTITIES ON INTEGER VACUUM MESH CSV :                       *
*                                                                     *
*           DG11LV = (D(VECT) / DS    . D(VECT) / DS  ) / JACOBIAN    *
*           DG22LV = (D(VECT) / DCHI  . D(VECT) / DCHI) / JACOBIAN    *
*           DG33LV = (D(VECT) / DPHI  . D(VECT) / DPHI) / JACOBIAN    *
*           DG12LV = (D(VECT) / DS    . D(VECT) / DCHI) / JACOBIAN    *
*                                                                     *
*       QUANTITIES ON HALF-INTEGER VACUUM MESH CSMV :                 *
*                                                                     *
*            ....M = ................AS ABOVE.....                    *
*                                                                     *
***********************************************************************
C
         COMMON /COMVAC/
     C   DG11LV,  DG22LV,   DG33LV,  DG12LV,
     C   DG11LMV, DG22LMV,  DG33LMV, DG12LMV,
     C   DGRMLV,  DGZMLV,   DGRMLMV, DGZMLMV
C
         COMPLEX
     C   DG11LV(NPV1,MPSMAX),    DG22LV(NPV1,MPSMAX),
     C   DG33LV(NPV1,MPSMAX),    DG12LV(NPV1,MPSMAX),
     C   DGRMLV(NPV1,MPSMAX),    DGZMLV(NPV1,MPSMAX),
     C   DG11LMV(NPV,MPSMAX),    DG22LMV(NPV,MPSMAX),
     C   DG33LMV(NPV,MPSMAX),    DG12LMV(NPV,MPSMAX),
     C   DGRMLMV(NPV,MPSMAX),    DGZMLMV(NPV,MPSMAX)




