C*COMDECK COMPLO
C
***********************************************************************
*                                                                     *
*    C.15  AUXILIARY RESULTS FOR PLOTS UNIQUELY                       *
*                                                                     *
*                                                                     *
*    RRCURV  : R POSITIONS OF INTERSECTIONS BETWEEN CONSTANT FLUX     *
*              LINE AND MAGNETIC FIELD ZERO CURVATURE LINE            *
*    RSHEAR  : LOCAL SHEAR ON (CHI, PSI) NODES OF ERATO MESH          *
*    RZCURV  : Z POSITIONS OF INTERSECTIONS BETWEEN CONSTANT FLUX     *
*              LINE AND MAGNETIC FIELD ZERO CURVATURE LINE            *
*                                                                     *
*    NCURV   : # OF CONSTANT FLUX SURFACES CROSSED BY THE MAGNETIC    *
*              FIELD ZERO CURVATURE LINE                              *
*                                                                     *
***********************************************************************
C
         COMMON /COMPLO/
     R     RRCURV  ,RSHEAR  ,RZCURV  ,
     I     NCURV
C
         DIMENSION
     R     RRCURV(4*NPPSI1),  RSHEAR(NPCHI,NPPSI1),  RZCURV(4*NPPSI1)

