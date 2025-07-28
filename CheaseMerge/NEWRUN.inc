C*COMDECK NEWRUN
C
***********************************************************************
*                                                                     *
*     C.21   ERATO NAMELIST                                           *
*                                                                     *
***********************************************************************
C
         DIMENSION
     R   ANGLE(16),   QTILDA(NPPSI1),  QS(NPPSI1),  WALL(10)
C
         LOGICAL   NLDIAG(1:20),
     L   NLEINQ,   NLGREN,        NLSYM
C
C---------------------------------------------------------------------
C
       NAMELIST /NEWRUN/
     R   AL0,      ANGLE,    ARROW,    BETA,     BETAP,    CPSRF,
     R   EPSCON,   EPSMAC,   QIAXE,    QSURF,    REXT,     RITOT,
     R   TSURF,    WNTORE,   WALL,     RINOR,    QCYL,     P0,
     I   NAL0AUTO,
     I   NER,      NEGP,     ITEST,    MEQ,      NCHI,     NDES,
     I   NFIG,     NITMAX,   NPRNT,    NPSI,     NSAVE,    NUA1,
     I   NUA2,     NUB1,     NUB2,     NUX,      NUPL,     NUSG,
     I   NV,       NUWA,     NVAC,     NVIT,     NWALL,
     L   NLDIAG,   NLEINQ,   NLGREN,   NLSYM
C
