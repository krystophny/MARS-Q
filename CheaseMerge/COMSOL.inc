C*COMDECK COMSOL
C
***********************************************************************
*                                                                     *
*     C.16         RESOLUTION PARAMETERS                              *
*                                                                     *
*     CEPS    : MINIMUM INTERMEDIATE VARIATION ACCEPTED FOR THE ERROR *
*               BETWEEN 2 ITERATIONS OF THE VARIATIONAL PROBLEM       *
*     CPI     : PI                                                    *
*     EPSLON  : FINAL VARIATION ACCEPTED FOR THE ERROR BETWEEN 2      *
*               ITERATIONS OF THE VARIATIONAL PROBLEM                 *
*     relax   : fraction of old solution in nonlinear iteration       *
*     RESIDU  : ERROR BETWEEN 2 ITERATIONS OF THE VARIATIONAL PROBLEM *
*     RESMAP  : ERROR BETWEEN 2 ITERATIONS OVER CURRENT PROFILE       *
*     SCALE   : SCALING FACTOR USED TO NORMALIZE TOTAL TOROIDAL       *
*               CURRENT                                               *
*     scexp   : for nsttp = 3 this is power to which computed scale   *
*               is raised before it is applied to the AT-coefficients *
*               normally 0 < scexp <= 1.                              *
*     SCHECK  : RESIDU OF A * PSI - B                                 *
*                                                                     *
*     NCON    : CONVERGENCE FLAG FOR VARIATIONAL PROBLEM              *
*     NINMAP  : MAXIMUM NUMBER OF ITERATIONS OVER MAPPING             *
*     NINSCA  : MAXIMUM NUMBER OF ITERATIONS OVER NON-LINEARITY       *
*                                                                     *
***********************************************************************
C
         COMMON /COMSOL/
     R       CEPS     ,CPI      ,EPSLON   ,relax    ,RESIDU   ,
     R       RESMAP   ,SCALE    ,SCALAC   ,scexp    ,SCHECK   ,
     I       NCON     ,NINMAP   ,NINSCA   

