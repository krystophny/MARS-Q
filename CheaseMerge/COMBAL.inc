C*COMDECK COMBAL
C
***********************************************************************
*                                                                     *
*     C.01         BALOONING + MERCIER STABILITY STORAGE              *
*                  EXPERIMENTAL PROFILES                              *
*                                                                     *
*     ABAL   : MATRIX OF BALOONING PROBLEM                            *
*     CFBAL  : COEFFICIENT WHICH LIMITS PRESSURE PROFILE (SEE EQ. (41)*
*              IN PUBLICATION)                                        *
*     CFNRESS: COEFFICIENT USED TO RESCALE FUNC PROFILE FOR NOPT.NE.1 *
*     CFNRESSO:INPUT CFNRESS VALUE, AS CFNRESS RESET TO 1. AFTER      *
*              INPUT PARAMETERS HAVE BEEN RESCALED                    *
*     CHI0   : INTEGRATION CONSTANT IN BALLOONING PROBLEM (SEE G IN   *
*              EQ. (17) OF PUBLICATION                                *
*     CPRESS : COEFFICIENT USED TO RESCALE BALLOONING OPTIMIZED       *
*              PRESSURE PROFILE (LAMBDA_2 IN PUBLICATION)             *
*     CPRESSO: INPUT CPRESS VALUE, AS CPRESS RESET TO 1. AFTER        *
*              INPUT PARAMETERS HAVE BEEN RESCALED                    *
*     D2RFUN : 2ND DERIVATIVES OF TT',I* OR I_PARALLEL  OBTAINED BY   *
*              CUBIC SPLINE INTERPOLATION OF TT',I* OR I_PARALLEL.    *
*              THIS ARRAY IS USED IF PROFILES ARE SPECIFIED WITH      *
*              EXPERIMENTAL DATA                                      *
*     D2RPRF : 2ND DERIVATIVES OF P' OBTAINED BY CUBIC SPLINE         *
*              INTERPOLATION OF P'. THIS ARRAY IS USED IF P' IS       *
*              SPECIFIED WITH EXPERIMENTAL DATA                       *
*     D2RPRM : 2ND DERIVATIVES OF P' OBTAINED BY CUBIC SPLINE         *
*              INTERPOLATION OF P'. THIS ARRAY IS USED FOR THE        *
*              BALLOONING OPTIMIZATION.                               *
*     FCSM   : S VALUES WHERE EXPERIMENTAL PROFILES ARE DEFINED       *
*     HMERCR : H OF GLASSER, GREENE AND JOHNSON FOR RESISTIVE         *
*              INTERCHANGE CRITERION (SEE EQ. (21) IN PUBLICATION)    *
*     PCS    : INTEGER S-MESH FOR BALLOONING OPTIMIZATION             *
*     PCSM   : HALF-INTEGER S-MESH FOR BALLOONING OPTIMIZATION        *
*     RFUN   : EXPERIMENTAL TT',I* OR I_PARALLEL PROFILE              *
*              (SEE SECTION 6.4.2 IN PUBLICATION)                     *
*     RJ1...RJ6 : INTEGRALS INVOLVED IN LOCAL INTERCHANGE CRITERIA    *
*              (SEE SECTION 4.1, EQ.(22) IN PUBLICATION               *
*     RPRF   : EXPERIMENTAL P' PROFILE                                *
*              (SEE SECTION 6.4.2 IN PUBLICATION)                     *
*     RPRM   : P' PROFILE USED FOR BALLOONING OPTIMIZATION.           *
*     SMERCI : MERCIER STABLITY TEST (SEE EQ. (19) IN PUBLICATION):   *
*                  SMERCI > 0 -----> STABLE                           *
*                  SMERCI < 0 -----> INSTABLE                         *
*     SMERCR : RESISTIVE INTERCHANGE CRITERION (SEE EQ. (20) IN       *
*              PUBLICATION):                                          *
*                  SMERCR > 0 -----> STABLE                           *
*                  SMERCR < 0 -----> INSTABLE                         *
*     XLAMB  : CONTAINS LAMBDA ARRAY IN THE BALLOONING OPTIMIZATION   *
*              ALGORITHM IN SECTION 5.5 OF THE PUBLICATION            *
*     XPPRDF : CONTAINS THE DELTA P' ARRAY IN THE BALLOONING          *
*              OPTIMIZATION ALGORITHM IN SECTION 5.5 OF THE           *
*              PUBLICATION                                            *
*     X2SRCH : CONTAINS THE DELTA LAMBDA ARRAY IN THE BALLOONING      *
*              OPTIMIZATION ALGORITHM IN SECTION 5.5 OF THE           *
*              PUBLICATION                                            *
*                                                                     *
*     XPPRMN, XPPRMX, XP0,XP1,XP2,XP3,XP4 ARE WORKING ARRAY USED      *
*     DURING THE BALLOONING OPTIMIZATION                              *
*                                                                     *
*     NCBAL  : BALLOONING STABILITY TEST : = 0 -----> STABLE          *
*                                          > 0 -----> INSTABLE        *
*     NCBLNS : FOR AN ASYMETRIC PLASMA SHAPE, SUM(NCBLNS(CHI-0)=NCBAL *
*                           NCBLNS(CHI-0) = 0 ---> STABLE             *
*                           NCBLNS(CHI-0) > 0 ---> INSTABLE           *
*     NPPF   : NUMBER OF EXPERIMENTAL DATA POINTS USED FOR THE        *
*              DEFINITION OF THE PROFILES                             *
*     NPPR   : NUMBER OF POINTS USED FOR THE P' PROFILE DURING THE    *
*              BALLOONING OPTIMIZATION                                *
*     NSRCH  : FLAG USED TO DIRECT BALLOONING OPTIMIZATION            *
*     NTURN  : NUMBER OF 2*PI INTERVALS USED FOR THE INTEGRATION OF   *
*              THE BALLOONING PROBLEM                                 *
*     N2BAL  : ARRAY USED TO INDICATE IF PRESSURE PROFILE HAS BEEN    *
*              LIMITED ACCORDING TO EQ. (41) IN THE PUBLICATION       *
*                           N2BAL = 0 ---> P' UNMODIFIED              *
*                           N2BAL = 1 ---> P' HAS BEEN LIMITED        *
*                                                                     *
*     NP0, NP1,NP2,NP3 ARE WORKING ARRAYS USED DURING THE             *
*     BALLOONING OPTIMIZATION                                         *
*                                                                     *
***********************************************************************
C
         COMMON /COMBAL/
     R   ABAL,    CFBAL,   CFNRESS,  CFNRESSO, CHI0,     CPRESS,
     R   CPRESSO, D2RFUN,  D2RPRM,   D2RPPF,   FCSM,     HMERCR,   PCS,      
     R   PCSM,    RFUN,    RJ1,      RJ2,      RJ3,
     R   RJ4,     RJ5,     RJ5P,     RJ6,      RPRM,     
     R   RPPF,    SMERCI,  SMERCR,   XLAMB,    XPPRMN,   
     R   XPPRMX,  XPPRDF,  XP0,      XP1,      XP2,      
     R   XP3,     XP4,     X2SRCH,
     I   NCBAL,   NCBLNS,  NPPF,     NPPR,     NP0,     
     I   NP1,     NP2,     NP3,      NP4,      NSRCH,    
     I   NTURN,   N2BAL
C
         DIMENSION
     R   ABAL(NPPSBAL,2,2*NPTURN*NPCHI+2), CHI0(NPBLC0),   
     R   D2RFUN(NPISO),                    D2RPRM(NPISO),
     R   D2RPPF(NPISO),                    FCSM(NPISO),               
     R   HMERCR(NPISO),                    PCS(NPISO),                    
     R   PCSM(NPISO),                      RFUN(NPISO), 
     R   RJ1(2*NPISO),                     RJ2(2*NPISO),
     R   RJ3(2*NPISO),                     RJ4(2*NPISO),
     R   RJ5(2*NPISO),                     RJ5P(2*NPISO),
     R   RJ6(2*NPISO),                     RPRM(NPISO),    
     R   RPPF(NPISO),                      SMERCI(NPISO),
     R   SMERCR(NPISO),                    XLAMB(NPISO),       
     R   XPPRMN(NPISO),                    XPPRMX(NPISO), 
     R   XPPRDF(NPISO),                    XP0(NPISO),    
     R   XP1(NPISO),                       XP2(NPISO),      
     R   XP3(NPISO),                       XP4(NPISO),    
     R   X2SRCH(NPISO),  
     I   NCBAL(NPISO),                     NCBLNS(NPISO,NPBLC0),
     I   NP0(NPISO),                       NP1(NPISO), 
     I   NP2(NPISO),                       NP3(NPISO), 
     I   NP4(NPISO),                       N2BAL(NPISO)

