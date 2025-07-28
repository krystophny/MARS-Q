*COMDECK INTEGC
C
C.....STATEMENT FUNCTIONS TO DO OVERLAP INTEGRALS FOR GENERALISED
C.....FEM METHOD
C
      REAL*8          ZNORM,PTRAP
      REAL*8          HP,HM,HC
      COMPLEX*16    FF,FF0,FF1,FFM,FFP,FGM,FGP,GG,GF
      COMPLEX*16    EQ,EQIP,EQIM,EQM,EQMP,EQMM
C
      FF(EQ,EQMM,EQMP) = ZNORM*(
     &  0.25*(1.-PTRAP)*(HM*EQMM + HP*EQMP) + 0.5*PTRAP*(HM + HP)*EQ )
C
      FF0(EQ,EQM) = ZNORM*(0.25*(1.-PTRAP)*HM*EQM + 0.5*PTRAP*HM*EQ )
C
      FF1(EQ,EQM) = ZNORM*(0.25*(1.-PTRAP)*HP*EQM + 0.5*PTRAP*HP*EQ )
C
      FFM(EQMM) = ZNORM * 0.25*(1.-PTRAP)*HM*EQMM
      FFP(EQMP) = ZNORM * 0.25*(1.-PTRAP)*HP*EQMP
C
      FGM(EQ,EQMM) = ZNORM * 0.5*HM*(PTRAP*EQ + (1.-PTRAP)*EQMM)
C
      FGP(EQ,EQMP) = ZNORM * 0.5*HP*(PTRAP*EQ + (1.-PTRAP)*EQMP)
C
      GF(EQ,EQM)   = ZNORM * 0.5*HC*(PTRAP*EQ + (1.-PTRAP)*EQM)
C
      GG(EQM,EQIM,EQIP) = ZNORM*
     &                HC*((1.-PTRAP)*EQM + 0.5*PTRAP*(EQIM+EQIP))
C
C----------------------------------------------------------------------




