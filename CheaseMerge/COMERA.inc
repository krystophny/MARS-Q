C*COMDECK COMERA
C
***********************************************************************
*                                                                     *
*     C.05         ERATO QUANTITIES                                   *
*                                                                     *
*     (CNR1,CNZ1) = (DR/DS,DZ/DS) ON (CHI,PSI) NODES NEEDED BY ERATO  *
*     (CNR2,CNZ2) = (DR/DPSI,DZ/DPSI) ON (CHI,PSI) NODES              *
*     (CR,CZ)     = (R,Z) OF (CHI,PSI) NODES                          *
*     EQ13 = EQ(13) NEEDED FOR ERATO DIAGNOSTICS                      *
*     EQ22 = EQ(22) NEEDED FOR ERATO DIAGNOSTICS                      *
*     EQ24 = EQ(24) NEEDED FOR ERATO DIAGNOSTICS                      *
*                                                                     *
*     EQ: CONTAINS QUANTITES NEEDED BY ERATO. (SEE [1], TABLE 2)      *
*     FOR EVERY PSI(I)-CHI(J) NODE, I=1,..,NPSI, J=1,..,NCHI          *
*                                                                     *
*     EQ( 1) = CS(I)                                                  *
*     EQ( 2) = CHI(J) WITH GENERAL JACOBIAN                           *
*     EQ( 3) = CS(I+1)                                                *
*     EQ( 4) = CHI(J+1) WITH GENERAL JACOBIAN                         *
*     EQ( 5) = CSM(I)                                                 *
*     EQ( 6) = CHIM(J) WITH GENERAL JACOBIAN                          *
*     EQ( 7) = O.          (CONTAINS MASS DENSITY)                    *
*     EQ( 8) = GAMMA * PRESSURE(I) / ( Q(MAGN. AXIS) * PSI-BOUND )    *
*     EQ( 9) = TOROIDAL FLUX T(I)                                     *
*     EQ(10) = 0.          (CONTAINS HELICAL PARAMETER)               *
*     EQ(11) = Q(I) / Q(MAGN. AXIS)                                   *
*     EQ(12) = PSI * R**2 / ( Q(MAGN. AXIS) * (GRAD-PSI)**2 )         *
*     EQ(13) = BETA-CHI-HAT DEFINED WITH GENERAL JACOBIAN             *
*     EQ(14) = R**2                                                   *
*     EQ(15) = 2 * ((D(LN(R**2)) / DS),CHI=CONST.)                    *
*     EQ(16) = ((D(LN(R**2)) / DCHI),PSI=CONST.)                      *
*     EQ(17) = H                                                      *
*     EQ(18) = K                                                      *
*     EQ(19) = R**2 / JACOBIAN                                        *
*     EQ(20) = ((D(EQ(19)) / DCHI),PSI=CONST.) / EQ(19)               *
*     EQ(21) = 2 * CSM(I) * PSI-BOUND *                               *
*              (D(INTEGRAL(JACOBIAN*TOROIDAL-FLUX*DCHI / R**2)) / DPSI)
*              WHERE THE INTEGRAL DOMAIN IS [O;CHI]                   *
*     EQ(22) = BETA-CHI-HAT WITH STRAIGHT FIELD LINES                 *
*     EQ(23) = 3 * NORM(GRAD-PSI) * (D(NORM(GRAD-PSI)) / DCHI)        *
*              WITH CONSTANT PSI                                      *
*     EQ(24) = CHI(J) WITH STRAIGHT FIELD LINES                       *
*     EQ(25) = CHI(J+1) WITH STRAIGHT FIELD LINES                     *
*                                                                     *
*     THE GENERAL JACOBIAN IS GIVEN BY                                *
*                                                                     *
*         J = CONST. * R**NER * NORM(GRAD-PSI)**NEGP                  *
*                                                                     *
*     THE STRAIGHT FIELD LINE JACOBIAN IS GIVEN BY                    *
*                                                                     *
*         J = CONST. * R**2                                           *
*                                                                     *
*     THE CONSTANT IS OBTAINED BY THE NORMALISATION OF THE CORRES-    *
*     PONDING CHI'S TO 2*PI. IT VARIES FROM ONE FLUX SURFACE TO       *
*     ANOTHER                                                         *
*                                                                     *
***********************************************************************
C
         COMMON /COMERA/
     R   CNR1   ,CNR2   ,CNZ1   ,CNZ2    ,CR     ,CZ     ,
     R   EQ13   ,EQ22   ,EQ24   ,EQ
C
         DIMENSION
     R   CNR1(NPCHI,NPPSI1),  CNR2(NPCHI,NPPSI1),   
     R   CNZ1(NPCHI,NPPSI1),  CNZ2(NPCHI,NPPSI1),   
     R   CR(NPCHI,NPPSI1),    CZ(NPCHI,NPPSI1),    
     R   EQ13(NPCHI,NPPSI1),  EQ22(NPCHI,NPPSI1),
     R   EQ24(NPCHI,NPPSI1),
     R   EQ(NDEQ,NPCHI,NPPSI1)
