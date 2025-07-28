C*COMDECK COMEQD
C**********************************************************************
C                                                                     *
C        COMMON FOR VALUES NEEDED FOR EQDSK TYPE OF OUTPUT (NIDEAL=6) *
C        P, P', G, GG' ON EQUIDISTANT PSI MESH 1..NRBOX               *
C        USE CPR, CPPR, TTP, TMF WITH NPSI1 SET TO NRBOX IN AUXVAL    *
C        PSI(I,J) ON (R,Z) GRID, EQUIDISTANT IN BOX:                  *
C        RBOXLFT TO RBOXLFT+RBOXLEN IN R (I=1,NRBOX)                  *
C        -ZBOXLEN/2. TO ZBOXLEN/2. IN Z  (J=1,NZBOX)                  *
C        THUS, ONE HAS TO EXTRAPOLATE PSI OUTSIDE PLASMA BOUNDARY     *
C                                                                     *
C        RBOXLFT, RBOXLEN, ZBOXLEN IN METERS                          *
C
C        R0EXP (METERS) AND B0EXP(TESLA) ARE USED TO RECONSTRUCT THE  *
C        EQUILIBRIUM IN MKSA VALUES. THEY EVOLVE IF NRSCAL=1, THUS ARE*
C        SAVED IN NOUT. THEY ARE READ FROM NIN IF NOPT=1, EXCEPT IF   *
C        THE VALUES ARE NEGATIVE IN NAMELIST.(SEE IODISK)             *
C                                                                     *
C**********************************************************************
C
         COMMON /COMEQD/
     R   EQDSPSI,  RBOXLFT, RBOXLEN, ZBOXLEN,
     I   NRBOX,    NZBOX
C
         DIMENSION
     R   EQDSPSI(NPISO,NPISO)
