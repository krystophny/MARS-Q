C*COMDECK COMMAP
C
***********************************************************************
*                                                                     *
*     C.12         MAPPING QUANTITIES                                 *
*                                                                     *
*     BCHIN    : BETA-CHI ON INTERSECTIONS BETWEEN FLUX SURFACES AND  *
*                CELL EDGES                                           *
*     BCHIO    : BETA-CHI WITH PRINCETON JACOBIAN ON INTERSECTIONS    *
*                BETWEEN FLUX SURFACES AND CELL EDGES                 *
*     CHIN     : CHI ON INTERSECTIONS BETWEEN FLUX SURFACES AND       *
*                CELL EDGES                                           *
*     CHIO     : CHI WITH PRINCETON JACOBIAN ON INTERSECTIONS         *
*                BETWEEN FLUX SURFACES AND CELL EDGES                 *
*     CHIOLD   : CHI WITH PRINCETON JACOBIAN ON ERATO CHI VALUES AT   *
*                PLASMA BOUNDARY                                      *
*     RHOVAC   : RHO VALUES OF (CHI,PSI) INTERSECTIONS ON PLASMA      *
*                SURFACE (FOR LION ONLY)                              *
*     RHOVACM  : RHO VALUES OF (CHIM,PSI) INTERSECTIONS ON PLASMA     *
*                SURFACE (FOR LION ONLY)                              *
*     SIGCHI   : SIGMA CORRESPONDING TO ERATO (CHI,CSM) NODES         *
*     SIGMAP   : SIGMA ON INTERSECTIONS BETWEEN FLUX SURFACES AND     *
*                CELL EDGES                                           *
*     SIGPEN   : EQUILIBRIUM SIGMA CORRESPONDING TO PENN MESH         *
*     TETCHI   : THETA CORRESPONDING TO ERATO (CHI,CSM) NODES         *
*     TETMAP   : THETA ON INTERSECTIONS BETWEEN FLUX SURFACES AND     *
*                CELL EDGES                                           *
*     TETPEN   : EQUILIBRIUM THETA CORRESPONDING TO PENN MESH         *
*     TETVAC   : THETA VALUES OF (CHI,PSI) INTERSECTIONS ON PLASMA    *
*                SURFACE (FOR LION ONLY)                              *
*     TETVACM  : THETA VALUES OF (CHIM,PSI) INTERSECTIONS ON PLASMA   *
*                SURFACE (FOR LION ONLY)                              *
*                                                                     *
***********************************************************************
C
         COMMON /COMMAP/
     R   BCHIN    ,BCHIO    ,CHIN     ,CHIO     ,
     R   CHIOLD   ,DRHOPI   ,RHOVAC   ,RHOVACI  ,
     R   RHOVACM  ,SIGCHI   ,SIGMAP   ,SIGPEN   ,
     R   TETCHI   ,TETMAP   ,TETPEN   ,TETVAC   ,
     R   TETVACI  ,TETVACM
C
         DIMENSION
     R   BCHIN(NTP2,2*NPISO),      BCHIO(NTP2,2*NPISO),
     R   CHIN(NTP2,2*NPISO),       CHIO(NTP2,2*NPISO),
     R   CHIOLD(NPCHI1),           DRHOPI(NPCHI,3),
     R   RHOVAC(NPCHI1),           RHOVACI(NPCHI,3),         
     R   RHOVACM(NPCHI1),          SIGCHI(NPCHI,NPPSI1),     
     R   SIGMAP(NTP2,2*NPISO),     SIGPEN(NPCHI*(NPMGS+1)),
     R   TETCHI(NPCHI,NPPSI1),     TETMAP(NTP2,2*NPISO),    
     R   TETPEN(NPCHI*(NPMGS+1)),  TETVAC(NPCHI1),           
     R   TETVACI(NPCHI,3),         TETVACM(NPCHI1)

