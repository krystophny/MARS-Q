C*COMDECK COMBND
C
***********************************************************************
*                                                                     *
*     C.03         BOUNDARY PARAMETERS                                *
*                                                                     *
*     ALZERO: NSURF=7:BOUNDARY AS SIN,COS: CENTRE IS (RC,RZ0C) AND    *
*     BPSCOS: RHO=ALZERO + SUM(BPSCOS(K)*COS(K*THET) +                *
*     BPSSIN:                  BPSSIN(K)*SIN(K*THET)), K=1,NFOURPB    *
*     BPS:   CONTAINS PARAMETERS TO SET UP PLASMA BOUNDARY SHAPE      *
*     BPSO:  AUXILIARY ARRAY USE TO STORE BPS                         *
*     RRBPS: R COORDINATES OF EXPERIMENTAL BOUNDARY POINTS            *
*     RZBPS: Z COORDINATES OF EXPERIMENTAL BOUNDARY POINTS            *
*     D2RBPS:2ND DERIVATIVE OF RRBPS WITH RESPECT TO THETA OBTAINED   *
*            BY CUBIC SPLINE INTERPOLATION                            *
*     D2ZBPS:2ND DERIVATIVE OF RZBPS WITH RESPECT TO THETA OBTAINED   *
*            BY CUBIC SPLINE INTERPOLATION                            *
*     NFOURPB: NUMBER OF FOURIER COMPONENTS USED TO SMOOTH            *
*              EXPERIMENTAL PLASMA BOUNDARY                           *
*     NBPS:  NUMBER OF EXPERIMENTAL BOUNDARY POINTS                   *
*                                                                     *
*     NBPSOUT, RRBPSOU, RZBPSOU: VALUES FOR OUTPUT ON EXPEQ.OUT       *
*                                AND EQDSK (COMPUT. IN RZBOUND)       *
*                                                                     *
***********************************************************************
C
         COMMON /COMBND/
     R   ALZERO, BPS,    BPSO,    BPSCOS, BPSSIN,
     R   D2RBPS, D2ZBPS, D2CBPS,  RRBPS,   RZBPS,  CNDRZ, RRBPSOU, 
     R   RRBPSNEW, RZBPSNEW,
     R   RZBPSOU,RZ0C,   TETBPS,  ZONE,   ZZERO, RW, RCOIL, HEND,
     I   NFOURPB,NBPS,   NBPSOUT, NI,     NF,
     I   NWBPS,  NW,     NDATA
C
         DIMENSION  BPS(14),       BPSO(14),
     R              BPSCOS(NPBPS),       
     R              BPSSIN(NPBPS),
     R              D2RBPS(NPBPS,NWBPS0), D2ZBPS(NPBPS,NWBPS0),
     R              D2CBPS(NPBPS,NWBPS0), 
     R              RRBPS(NPBPS,NWBPS0),  RZBPS(NPBPS,NWBPS0),
     R              RRBPSNEW(NPBPS),  RZBPSNEW(NPBPS),
     R              CNDRZ(NPBPS,NWBPS0),
     R              RRBPSOU(NPBPS),RZBPSOU(NPBPS),
     R              TETBPS(NPBPS,NWBPS0), RW(NWBPS0), RCOIL(5),
     I              NW(NWBPS0)


