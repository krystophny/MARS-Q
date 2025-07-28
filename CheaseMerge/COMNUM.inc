C*COMDECK COMNUM
C
***********************************************************************
*                                                                     *
*     C.13         MESH NUMBERS                                       *
*                                                                     *
*     MSMAX  : # OF POSITIF POLOIDAL MODES FOR LINEAR RESISTIF        *
*              STABILITY CODE                                         *
*     NBAND  : BAND WIDTH OF MATRIX A                                 *
*     NBLC0  : # OF INTEGRATION CONSTANTS USED FOR BALOONING          *
*              CALCULATION                                            *
*     NCHI   : # OF CHI INTERVALS                                     *
*     NCHI1  : NCHI + 1                                               *
*     NDIM   : CONTAINS FINAL NS,NT AND NISO VALUES DURING ITERATIONS *
*              OVER MESH                                              *
*     NISO   : # OF INTERVALS OF S-IPRIME MESH                        *
*     NISOO  : OLD VALUE OF NISO DURING ITERATION OVER MESH           *
*     NMGAUS : # OF INTEGRATION POINTS IN THE THETA DIRECTION FOR THE *
*              LINE INTEGRALS COMPUTED IN THE MAPPINGS                *
*     NPSI   : # OF S INTERVALS                                       *
*     NPSI1  : NPSI + 1                                               *
*     NS     : # OF SIGMA INTERVALS                                   *
*     NSGAUS : # OF INTEGRATION POINTS IN THE SIGMA DIRECTION         *
*     NSMAX  : # OF POSITIF TOROIDAL MODES FOR LINEAR RESISTIF        *
*              STABILITY CODE                                         *
*     NSTMAX : NS1 * NT                                               *
*     NSO    : OLD VALUE OF NS DURING ITERATION OVER MESH             *
*     NS1    : NS + 1                                                 *
*     NT     : # OF THETA INTERVALS                                   *
*     NTGAUS : # OF INTEGRATION POINTS IN THE THETA DIRECTION         *
*     NTO    : OLD VALUE OF NT DURING ITERATION OVER MESH             *
*     NT1    : NT + 1                                                 *
*     NT2    : NT + 2                                                 *
*     NV     : # OF RADIAL INTERVALS FOR VACUUM DISCRETIZATION (MARS) *
*     NV1    : NV + 1                                                 *
*     NWGAUS : NSGAUS * NTGAUS                                        *
*     N4NSNT : 4 * NSTMAX                                             *
*                                                                     *
***********************************************************************
C
         COMMON /COMNUM/
     I   MSMAX,   NBAND,   NBLC0,   NCHI,    NCHI1,   NDIM,    
     I   NISO,    NISOO,   NMGAUS,  NPSI,    NPSI1,   NS,      
     I   NSGAUS,  NSMAX,   NSTMAX,  NSO,     NS1,     NT,      
     I   NTGAUS,  NTNOVA,  NTO,     NT1,     NT2,     NV,     
     I   NV1,     NWGAUS,  N4NSNT
C
         DIMENSION
     I   NDIM(5)

