      MODULE gijlm
C
      COMPLEX*16,DIMENSION(:,:),ALLOCATABLE::
     C     DG11L, DG22L, DG33L, DG12L, GCHDZ, GSDZ, GBZ, GBR,
     C     DG11LM,DG22LM,DG33LM,DG12LM,GCHDZM,GSDZM,GBZM,GBRM
C-----------------------------------------------------------------------
C
C  METRIC TENSOR QUANTITIES:  (FUNCTIONS OF S AND CHI)
C
C            G11L  =  S (1)- S (1) TENSOR ELEMENT ON INTEGER MESH (CS)
C                                  L MEANS LOW INDICES
C            G22L  = CHI(2)-CHI(2) TENSOR ELEMENT ON INTEGER MESH (CS)
C            G33L  = PHI(3)-PHI(3) TENSOR ELEMENT ON INTEGER MESH (CS)
C            G12L  =  S (1)-CHI(2) TENSOR ELEMENT ON INTEGER MESH (CS)
C
C            ....M = ................AS ABOVE.....ON   HALF  MESH (CSM)
C
C           D..... = D*G...    WHERE D IS THE INVERSE OF THE JACOBIAN
C
C-----------------------------------------------------------------------
      END MODULE gijlm
