C*COMDECK COMVEV
C
***********************************************************************
*                                                                     *
*   C.19    VARIABLES FOR SOLOVEV CASE ONLY                           *
*                                                                     *
*                                                                     *
*     CPSI1T : THEORETICAL VALUES OF PSI                              *
*     DIFFDS : DIFFERENCE BETWEEN NUMERIC AND THEORETICAL DPSI / DS   *
*     DIFFDT : DIFFERENCE BETWEEN NUMERIC AND THEORETICAL DPSI / DT   *
*     DIFFP  : DIFFERENCE BETWEEN NUMERIC AND THEORETICAL PSI         *
*     DPDSNU : NUMERIC VALUES OF DPSI / DS                            *
*     DPDTNU : NUMERIC VALUES OF DPSI / DT                            *
*     DPDSTH : THEORETICAL VALUES OF DPSI / DS                        *
*     DPDTTH : THEORETICAL VALUES OF DPSI / DT                        *
*     RESDPS : RESIDU OF DPSI / DS                                    *
*     RESDPT : RESIDU OF DPSI / DT                                    *
*     RESDST : RESIDU OF D2PSI / (DS*DT)                              *
*     RESPSI : RESIDU OF PSI                                          *
*                                                                     *
***********************************************************************
C
         COMMON /COMVEV/
     R   CPSI1T  ,DIFFDS  ,DIFFDT  ,DIFFP   ,DIFFST  ,
     R   DPDSNU  ,DPDTNU  ,DPDSTH  ,DPDTTH  ,D2PSTN  ,
     R   D2PSTT  ,RESDPS  ,RESDPT  ,RESDST  ,RESPSI
C
         DIMENSION
     R   CPSI1T(NSNT),    DIFFDS(NSNT),    DIFFDT(NSNT),
     R   DIFFP(NSNT),     DIFFST(NSNT),    DPDSNU(NSNT),
     R   DPDTNU(NSNT),    DPDSTH(NSNT),    DPDTTH(NSNT),
     R   D2PSTN(NSNT),    D2PSTT(NSNT)

