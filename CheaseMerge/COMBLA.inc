C*COMDECK COMBLA
C
***********************************************************************
*                                                                     *
*     C.02         SOLUTION + MATRIX STORAGE                          *
*                                                                     *
*                                                                     *
*     A      : MATRIX OF VARIATIONNAL PROBLEM (SEE LEFT HAND SIDE     *
*              OF EQ. (27) IN PUBLICATION                             *
*     B      : RHS OF VARIATIONAL PROBLEM (SEE RIGHT HAND SIDE OF     *
*              EQ. (27) IN PUBLICATION                                *
*     CPSI   : SOLUTION VECTOR OF VARIATIONNAL PROBLEM. THE ORGANI-   *
*              SATION OF THE VARIABLES IS THE FOLLOWING               *
*                  1) PSI(I)                                          *
*                  2) DPSI / DS (I)                                   *
*                  3) DPSI / DT (I)                                   *
*                  4) D2PSI / DSDT (I)                                *
*              FOR I = 1,..,NSNT. THE NUMEROTATION IS ALTERNATE       *
*              UP-DOWN FOR THE MESH                                   *
*     CPSI0  : AUXILIARY ARRAY USED TO STORE CPSI                     *
*     CPSICL : SAME AS CPSI, BUT INVERSE CLOCKWISE NUMBERING          *
*     CPSILI : AS CPSI, BUT CONTAINS PSI_C, PSI_R AND PSI_Z IN        *
*              EQ. (30) OF PUBLICATION INSTEAD OF PSI, FIRST AND      *
*              CROSS DERIVATIVES OF PSI WITH RESPECT TO SIGMA AND     *
*              THETA COMPUTED WITH PSI_C, PSI_R AND PSI_Z             *
*                                                                     *
***********************************************************************
C
         COMMON /COMBLA/
     R   A       ,B       ,CPSI    ,CPSIO   ,CPSICL  ,CPSILI
C
         DIMENSION
     R   A(NPBAND,NP4NST),   B(NP4NST),         CPSI(NP4NST),
     R   CPSIO(NP4NST),      CPSICL(NP4NST),    CPSILI(NP4NST)


