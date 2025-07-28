C*COMDECK COMINT
C
***********************************************************************
*                                                                     *
*     C.08         INTEGRATION PARAMETERS                             *
*                                                                     *
*     CW    : 4-POINT GAUSS QUADRATURE WEIGHTS IN [0; 1] X [0; 1]     *
*     DZETA : LOCATION OF INTEGRATION POINTS IN [0; 1] X [0; 1]       *
*     FB    : VALUES OF BASIS FUNCTIONS ON INTEGRATION POINTS IN      *
*             (SIGMA; THETA) MESH                                     *
*     RSINT : SIGMA VALUES OF INTEGRATION POINTS                      *
*     RTINT : THETA VALUES OF INTEGRATION POINTS                      *
*     YDRSDT: D(YRST(J)) / DT ,J=1,...,16*NT                          *
*     YRST  : BOUNDARY VALUES ON RTINT VALUES                         *
*                                                                     *
*     MPLA  : HORIZONTAL LOCATIONS IN MATRIX A                        *
*     NPLAC : VERTICAL LOCATIONS IN MATRIX A                          *
*     NUPDWN: UP-DOWN NUMEROTAION AS A FUNCTION OF INVERSE CLOCKWISE  *
*             NUMEROTATION                                            *
*                                                                     *
***********************************************************************
C
         COMMON /COMINT/
     I   MPLA     ,NPLAC    ,NUPDWN   ,
     R   CW       ,DZETA    ,FB       ,RSINT    ,
     R   RTINT    ,YDRSDT   ,YRST
C
         DIMENSION
     I   MPLA(NPT,4,4),     NPLAC(NPT,4),      NUPDWN(NSNT),
     R   CW(NPGAUS),        DZETA(NPGAUS,2),   FB(NPT,NPGAUS,16,NPS),
     R   RSINT(NPS,NPGAUS), RTINT(NPT,NPGAUS), YDRSDT(NPT,NPGAUS),
     R   YRST(NPT,NPGAUS)
