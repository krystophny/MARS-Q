C*COMDECK COMESH
C
***********************************************************************
*                                                                     *
*     C.06         MESH QUANTITIES                                    *
*                                                                     *
*     NMAG     : CELL CONTAINING MAGNETIC AXIS                        *
*     NMESHA   : 0 -----> NO PACKING FOR S-MESH                       *
*                1 -----> PACKING FOR S-MESH ON SPECIFIED S SURFACES  *
*                2 -----> PACKING FOR S-MESH ON SPECIFIED Q SURFACES  *
*     NMESHB   : 0 -----> NO PACKING FOR S-IPRIME MESH                *
*                1 -----> PACKING FOR S-IPRIME MESH                   *
*     NMESHC   : 0 -----> NO PACKING FOR SIGMA-MESH                   *
*                1 -----> PACKING FOR SIGMA-MESH                      *
*     NMESHD   : 0 -----> NO PACKING FOR THETA-MESH                   *
*                1 -----> PACKING FOR THETA-MESH                      *
*     NMESHE   : 0 -----> NO PACKING FOR CHI-MESH                     *
*                1 -----> PACKING FOR CHI-MESH                        *
*     NPOIDA   : # OF LORENTZIANS USED FOR S-MESH PACKING             *
*     NPOIDB   : # OF LORENTZIANS USED FOR S-IPRIME MESH PACKING      *
*     NPOIDC   : # OF LORENTZIANS USED FOR SIGMA-MESH PACKING         *
*     NPOIDD   : # OF LORENTZIANS USED FOR THETA-MESH PACKING         *
*     NPOIDE   : # OF LORENTZIANS USED FOR CHI-MESH PACKING           *
*     NPOIDQ   : # OF RATIONNAL Q SURFACES WHERE THE S-MESH IS PACKED *
*                                                                     *
*     APLACE   : LOCATIONS OF LORENZTIANS FOR S-MESH PACKING          *
*     AWIDTH   : WIDTH OF LORENTZIANS FOR S-MESH PACKING              *
*     BPLACE   : LOCATIONS OF LORENZTIANS FOR S-IPRIME MESH PACKING   *
*     BWIDTH   : WIDTH OF LORENTZIANS FOR S-IPRIME MESH PACKING       *
*     B0EXP    : EXPERIMENTAL VALUE OF B0 IN TESLA (FOR EXPEQ.OUT)    *
*     CHI      : INTEGER CHI MESH                                     *
*     CHIM     : HALF INTEGER CHI MESH                                *
*     CPLACE   : LOCATIONS OF LORENZTIANS FOR SIGMA MESH PACKING      *
*     CS       : INTEGER S MESH                                       *
*     CSIG     : SIGMA MESH                                           *
*     CSIGO    : STORAGE SPACE FOR OLD SIGMA MESH DURING ITERATION    *
*                OVER MESH                                            *
*     CSIPR    : S-MESH FOR I-PRIME DEFINITION (HALF-INTEGER)         *
*     CSIPRI   : S-MESH FOR I-PRIME DEFINITION (INTEGER)              *
*     CSIPRO   : STORAGE SPACE FOR OLD CSIPR MESH DURING ITERATION    *
*                OVER MESH                                            *
*     CSM      : HALF-INTEGER S MESH                                  *
*     CSMV     : HALF-INTEGER VACUUM S MESH FOR MARS                  *
*     CSPEN    : S MESH FOR PENN                                      *
*     CSV      : INTEGER VACUUM S MESH FOR MARS                       *
*     CT       : THETA MESH                                           *
*     CTO      : STORAGE SPACE FOR OLD THETA MESH DURING ITERATION    *
*                OVER MESH                                            *
*     CTPEN    : THETA MESH FOR PENN                                  *
*     CTXT     : THETA MESH FOR XTOR                                  *
*     CWIDTH   : WIDTH OF LORENTZIANS FOR SIGMA MESH PACKING          *
*     DPLACE   : LOCATIONS OF LORENZTIANS FOR THETA MESH PACKING      *
*     DRSDT    : STORAGE SPACE FOR D(RHOS(CT(J))) / DT , J=1,...,NT   *
*     DWIDTH   : WIDTH OF LORENTZIANS FOR THETA MESH PACKING          *
*     EPLACE   : LOCATIONS OF LORENZTIANS FOR CHI MESH PACKING        *
*     EWIDTH   : WIDTH OF LORENTZIANS FOR CHI MESH PACKING            *
*     PSISCL   : FRACTION OF POLOIDAL FLUX KEPT FOR STABILITY         *
*                CALCULATION                                          *
*     QPLACE   : Q VALUES WHERE THE S MESH IS PACKED                  *
*     QWIDTH   : WIDTH OF LORENTZIANS FOR S MESH PACKING FOR NMESHA=2 *
*                OPTION                                               *
*     R0       : R OF MESH AXIS                                       *
*     R0EXP    : EXPERIMENTAL VALUE OF R0 IN METERS (FOR EXPEQ.OUT)   *
*     R0O      : STORE R0                                             *
*     R0W      : R OF VACUUM MESH AXIS FOR MARS                       *
*     R0WO     : STORE R0W                                            *
*     RC       : R REFERENCE USED TO COMPUTE PLASMA BOUNDARY          *
*     REXT     : DISTANCE BETWEEN PLASMA AND WALL                     *
*     RHOS     : PLASMA BOUNDARY OF THETA = CT(J), J=1,...,NT         *
*     RZ0      : Z OF MESH AXIS                                       *
*     RZ0O     : STORE RZ0                                            *
*     RZ0W     : Z OF VACUUM MESH AXIS FOR MARS                       *
*     RZ0WO    : STORE RZ0W                                           *
*     SOLPDA   : FRACTION OF UNDENSIFIED S-MESH                       *
*     SOLPDB   : SAME AS SOLPDA, BUT FOR CSIPR-MESH                   *
*     SOLPDC   : SAME AS SOLPDA, BUT FOR SIGMA-MESH                   *
*     SOLPDD   : SAME AS SOLPDA, BUT FOR THETA-MESH                   *
*     SOLPDE   : SAME AS SOLPDA, BUT FOR CHI-MESH                     *
*                                                                     *
***********************************************************************
C
         COMMON /COMESH/
     R   APLACE   ,AWIDTH   ,BPLACE   ,BWIDTH   ,B0EXP    ,
     R   CHI      ,CHIM     ,CPLACE   ,CS       ,CSIG     ,
     R   CSIGO    ,CSIPR    ,CSIPRI   ,CSIPRO   ,CSM      ,
     R   CSMV     ,CSPEN    ,CSV      ,CT       ,CTO      ,
     R   CTPEN    ,CWIDTH   ,CTXT     ,DPLACE   ,DRSDT    ,
     R   DWIDTH   ,EPLACE   ,EWIDTH   ,PSISCL   ,QPLACE   ,
     R   QWIDTH   ,R0       ,R0EXP    ,R0O      ,R0W      ,
     R   R0WO     ,RC       ,REXT     ,RHOS     ,RZ0      ,
     R   RZ0O     ,RZ0W     ,RZ0WO    ,SOLPDA   ,SOLPDB   ,
     R   SOLPDC   ,SOLPDD   ,SOLPDE   ,ROTE     ,QWIDTH0  ,
     R   VZ2GP    ,ZDEL     ,ZFRC     ,ZPTS,
     R   R0EXP0   ,B0EXP0,
     I   NMAG     ,NMESHA   ,NMESHB   ,NMESHC   ,NMESHD   ,
     I   NMESHE   ,NPOIDA   ,NPOIDB   ,NPOIDC   ,NPOIDD   ,
     I   NPOIDE   ,NPOIDQ   ,NTOR     ,NVACUUMRNW, KMETHOD,
     I   NPTS
C
         DIMENSION
     R   APLACE(NPPACK),         AWIDTH(NPPACK), BPLACE(NPPACK),
     R   BWIDTH(NPPACK),         CHI(NPCHI1),    CHIM(NPCHI),
     R   CPLACE(NPPACK),         CS(NPPSI1),     CSIG(NSP1),
     R   CSIGO(NSP1),            CSIPR(NPISO),   CSIPRI(NPISO),
     R   CSIPRO(NPISO),          CSM(NPPSI1),    CSMV(NPV),      
     R   CSPEN(NPPSI1*(NPMGS+1)),CSV(NPV1),      CT(NTP1),       
     R   CTPEN(NPCHI*(NPMGS+1)), CTO(NTP1),      CWIDTH(NPPACK),     
     R   CTXT(NPCHI1),           DPLACE(NPPACK), DRSDT(NTP1),    
     R   DWIDTH(NPPACK),         EPLACE(NPPACK), EWIDTH(NPPACK),     
     R   QPLACE(NPPACK),         QWIDTH(NPPACK), RHOS(NTP1),
     R   ZPTS(10)
