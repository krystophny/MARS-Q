C*COMDECK COMISO
C
***********************************************************************
*                                                                     *
*     C.10         QUANTITIES ON GAUSSIAN QUADRATURE POINTS ON        *
*                  CONSTANT FLUX SURFACES                             *
*                                                                     *
*                                                                     *
*     BCHISO  : BETA-CHI ON GAUSS QUADRATURE POINTS ALONG CONSTANT    *
*               FLUX SURFACES                                         *
*     BNDISO  : PLASMA BOUNDARY ON GAUSS QUADRATURE POINTS ALONG      *
*               CONSTANT FLUX SURFACES                                *
*     CHIISO  : CHI ON GAUSS QUADRATURE POINTS ALONG CONSTANT         *
*               FLUX SURFACES                                         *
*     DGNISO  : NORMAL DERIVATIVE OF (GRAD-PSI)**2  ON GAUSS QUADRA-  *
*               TURE POINTS ALONG CONSTANT FLUX SURFACES              *
*     DPSISO  : D(PSI) / D(SIGMA) ON GAUSS QUADRATURE POINTS ALONG    *
*               CONSTANT FLUX SURFACES                                *
*     DPTISO  : D(PSI) / D(THETA) ON GAUSS QUADRATURE POINTS ALONG    *
*               CONSTANT FLUX SURFACES                                *
*     DRNISO  : NORMAL DERIVATIVE OF R ON GAUSS QUADRATURE POINTS     *
*               ALONG CONSTANT FLUX SURFACES                          *
*     GPISO   : (GRAD-PSI)**2 ON GAUSS QUADRATURE POINTS ALONG        *
*               CONSTANT FLUX SURFACES                                *
*     RHOISO  : RHO ON GAUSS QUADRATURE POINTS ALONG CONSTANT FLUX    *
*               SURFACES                                              *
*     RRISO   : R ON GAUSS QUADRATURE POINTS ALONG CONSTANT FLUX      *
*               SURFACES                                              *
*     RZISO   : Z ON GAUSS QUADRATURE POINTS ALONG CONSTANT FLUX      *
*               SURFACES                                              *
*     SIGPSI  : SIGMA ON GAUSS QUADRATURE POINTS ALONG CONSTANT FLUX  *
*               SURFACES                                              *
*     TETPSI  : THETA ON GAUSS QUADRATURE POINTS ALONG CONSTANT FLUX  *
*               SURFACES                                              *
*     WGTPSI  : INTEGRATION WEIGHTS ON GAUSS QUADRATURE POINTS        *
*               ALONG CONSTANT FLUX SURFACES                          *
*                                                                     *
***********************************************************************
C
         COMMON /COMISO/
     R   BCHISO  ,BNDISO  ,CHIISO  ,DGNISO  ,DGRISO  ,
     R   DGZISO  ,DPSISO  ,DPTISO  ,DPRISO  ,DPZISO  ,
     R   DRNISO  ,GPISO   ,RHOISO  ,RRISO   ,RZISO   ,
     R   SIGPSI  ,TETPSI  ,WGTPSI
C
         DIMENSION
     R   BCHISO(NPMGS*NTP1),         BNDISO(NPMGS*NTP1,2*NPISO),
     R   CHIISO(NPMGS*NTP1),         DGNISO(NPMGS*NTP1,2*NPISO),
     R   DGRISO(NPMGS*NTP1,2*NPISO), DGZISO(NPMGS*NTP1,2*NPISO),
     R   DPSISO(NPMGS*NTP1,2*NPISO), DPTISO(NPMGS*NTP1,2*NPISO),
     R   DPRISO(NPMGS*NTP1,2*NPISO), DPZISO(NPMGS*NTP1,2*NPISO),
     R   DRNISO(NPMGS*NTP1,2*NPISO), GPISO(NPMGS*NTP1,2*NPISO),
     R   RHOISO(NPMGS*NTP1,2*NPISO), RRISO(NPMGS*NTP1,2*NPISO),
     R   RZISO(NPMGS*NTP1,2*NPISO),  SIGPSI(NPMGS*NTP1,2*NPISO),
     R   TETPSI(NPMGS*NTP1,2*NPISO), WGTPSI(NPMGS*NTP1,2*NPISO)

