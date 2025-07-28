C*COMDECK COMPHY
C
***********************************************************************
*                                                                     *
*     C.14         PHYSICAL PARAMETERS                                *
*                                                                     *
*     AFBS     : PRESCRIPTION OF BOOTSTRAP CURRENT (NPP = 1)          *
*     AFBS2    : PRESCRIPTION OF BOOTSTRAP CURRENT (NPP = 2)          *
*     AP       : PARAMETERS PRESCRIBING P-PRIME(PSI) PROFILE (NPP = 1)*
*     AP2      : PARAMETERS PRESCRIBING P-PRIME(S)   PROFILE (NPP = 2)*
*     AREA     : SURFACE OF THE PLASMA SECTION                        *
*     ASPCT    : INVERSE ASPECT RATIO                                 *
*     AT       : PRESCRIPTION OF TT', I* AND I_PARALLEL               *
*     AT2      : PRESCRIPTION OF TT', I* AND I_PARALLEL (NIPR=2)      *
*     AT3      : PRESCRIPTION OF TT', I* AND I_PARALLEL (NIPR=3)      *
*     AT4      : PRESCRIPTION OF TT', I* AND I_PARALLEL (NIPR=4)      *
*     BEANS    : INDENTATION OF PLASMA SECTION                        *
*     BETA     : TOROIDAL BETA OF EQUILIBRIUM                         *
*     BETAP    : POLOIDAL BETA OF EQUILIBRIUM                         *
*     BETAS    : BETA-STAR OF EQUILIBRIUM                             *
*     BETAX    : EXPERIMENTAL TOROIDAL BETA OF EQUILIBRIUM            *
*     BSFRAC   : FRACTION OF BOOTSTRAP CURRENT                        *
*     CETA     : PLASMA SURFACE DEFINITION                            *
*     CONVF    : CONVERGENCE FACTOR BETWEEN DIFFERENT DEFINITIONS     *
*                OF BETA, BETA_POL, LI                                *
*     CPBAR    : AVERAGED PRESSURE OF EQUILIBRIUM                     *
*     CPDP0    : D(CP) / D(PSI) ON MAGNETIC AXIS (CP IS THE CONSTANT  *
*                INVOLVED IN THE JACOBIAN WHICH DEFINES THE CHI-MESH) *
*     CPP      : P-PRIME FOR SOLOVEV CASE                             *
*     CPPF     : PRESSURE PEAKING FACTOR                              *
*     CPSI0    : CP ON MAGNETIC AXIS (CP IS THE CONSTANT INVOLVED IN  *
*                THE JACOBIAN WHICH DEFINES THE CHI-MESH)             *
*     CPSRF    : PSI AT PLASMA EDGE AFTER RENORMALIZATION OF EQUILI-  *
*                BRIUM FOR STABILITY CODES                            *
*     CQ0      : Q ON MAGNETIC AXIS USED TO DEFINE A SOLOVEV CASE     *
*     CSSPEC   : S VALUE WHERE Q HAS TO BE SPECIFIED                  *
*     CUROLD   : TOTAL TOROIDAL CURRENT OF EQUILIBRIUM BEFORE         *
*                RENORMALIZATIONS                                     *
*     CURRT    : TOTAL TOROIDAL CURRENT IMPOSED TO AVOID TRIVIAL      *
*                SOLUTION OF VARIATIONAL PROBLEM                      *
*     DELTA    : DELTA FOR X-POINT PLASMA SECTION                     *
*     DENS0    : PLASMA DENSITY ON MAGNETIC AXIS                      *
*     DPDP0    : P-PRIME ON MAGNETIC AXIS                             *
*     DQDP0    : Q-PRIME ON MAGNETIC AXIS                             *
*     DTTP0    : T*T-PRIME ON MAGNETIC AXIS                           *
*     ELONG    : ELONGATION                                           *
*     ETAI     : D(LOG T_E)/D(LOG N). SEE HIRSHMAN DEFINITION OF      *
*                BOOTSTRAP CURRENT                                    *
*     GAMMA    : RATIO OF SPECIFIC HEATH                              *
*     CP0      : PRESSURE ON AXIS                                     *
*     PANGLE   : ANGLE FOR PROFILE DEFINITION IN CARTESIAN PLANE      *
*     PREDGE   : PRESSURE AT PLASMA EDGE                              *
*     PSI0     : MINIMUM OF PSI OF EQUILIBRIUM ( < 0 )                *
*     QCYL     : CYLINDRICAL Q                                        *
*     QSPEC    : VALUE OF Q SPECIFIED AT CSSPEC                       *
*     Q0       : Q AT MAGNETIC AXIS                                   *
*     RFPF     : F FOR REVERSED PINCH                                 *
*     RFPT     : THETA FOR REVERSED PINCH                             *
*     RIBSNOR  : BOOTSTRAP CURRENT WITH ERATO NORMALIZATION           *
*     RINOR    : TOTAL CURRENT WITH ERATO NORMALIZATION               *
*     RIPR0    : I-PRIME ON AXIS                                      *
*     RITBS    : TOTAL BOOTSTRAP CURRENT WITHOUT COLLISIONALITY       *
*     RITBSC   : TOTAL BOOTSTRAP CURRENT INCLUDING COLLISIONALITY     *
*     RITOT    : TOTAL CURRENT AFTER SCALE OF EQUILIBRIUM BY NOREPT   *
*     RJDTB0   : <J.B> AT MAGNETIC AXIS                               *
*     RMAG     : R OF MAGNETIC AXIS                                   *
*     RMAGO    : STORE RMAG                                           *
*     RNU      : NU EXPONENT FOR X-POINT PLASMA SECTION               *
*     RRAXIS   : R OF MAGNETIC AXIS                                   *
*     RZAXIS   : Z OF MAGNETIC AXIS                                   *
*     RZION    : ION CHARGE IN HRSHAMN DEFINITION OF BOOTSTRAP        *
*                CURRENT                                              *
*     RZMAG    : Z OF MAGNETIC AXIS                                   *
*     RZMAGO   : STORE RZMAG                                          *
*     SCALNE   : NPROFZ = 1 -----> DENSITY = SCALNE * NORMALIZED      *
*                                  PROFILE                            *
*     SGMA     : PLASMA BOUNDARY DEFINITION                           *
*     SMAG     : SIGMA OF MAGNETIC AXIS                               *
*     SPSI0    : ABSOLUTE VALUE OF MINIMUM OF PSI FOR SOLOVEV CASE    *
*     SPSIM    : MINIMUM OF PSI : SPSIM < 0 FOR EQUILIBRIUM           *
*                                 SPSIM = 0 FOR ERATO                 *
*     SPSIMO   : STORE SPSIM                                          *
*     TEMP0    : PLASMA TEMPERATURE ON MAGNETIC AXIS                  *
*     THEAT0   : PLASMA BOUNDARY DEFINITION                           *
*     TMAG     : THETA OF MAGNETIC AXIS                               *
*     TRIANG   : PLASMA TRIANGULARITY                                 *
*     TRIPLT   : PLASMA BOUNDARY DEFINITION                           *
*     T0       : TOROIDAL FLUX ON MAGNETIC AXIS                       *
*     VOLUME   : PLASMA VOLUME                                        *
*     XI       : PLASMA BOUNDARY DEFINITION                           *
*     WMAGP    : POLOIDAL MAGNETIC FIELD ENERGY                       *
*                                                                     *
*     NSOUR    : MAXIMUM DEGREE OF POLYNOMS USED TO DEFINE P-PRIME AND*
*                T-TPRIME PROFILES FOR NPP = NSTTP = 0                *
*                                                                     *
***********************************************************************
C
         COMMON /COMPHY/
     R   AFBS     ,AFBS2    ,AP       ,AP2      ,AREA     ,
     R   ASPCT    ,AT       ,AT2      ,AT3      ,AT4      ,
     R   BEANS    ,BETA     ,BETAP    ,BETAS    ,BETAX    ,
     R   BSFRAC   ,CETA     ,CONVF    ,CPBAR    ,CPDP0    ,
     R   CPP      ,CPPF     ,CPSI0    ,CPSRF    ,CQ0      ,
     R   CSSPEC   ,CUROLD   ,CURRT    ,DELTA    ,DENS0    ,
     R   DPDP0    ,DQDP0    ,DTTP0    ,ELONG    ,ETAEI    ,
     R   GAMMA    ,CP0      ,PANGLE   ,PREDGE   ,PSI0     ,
     R   QCYL     ,QSPEC    ,Q0       ,RFPF     ,RFPT     ,
     R   RIBSNOR  ,RINOR    ,RIPR0    ,RITBS    ,RITBSC   ,
     R   RITOT    ,RJDTB0   ,RMAG     ,RMAGO    ,RNU      ,
     R   RRAXIS   ,RZAXIS   ,RZION    ,RZMAG    ,RZMAGO   ,
     R   SCALNE   ,SGMA     ,SMAG     ,SPSI0    ,SPSIM    ,
     R   SPSIMO   ,TEMP0    ,THETA0   ,TMAG     ,TRIANG   ,
     R   TRIPLT   ,T0       ,VOLUME   ,XI       ,WMAGP    ,
     R   QSHAVE   ,PSISHAVE ,   
     I   NSOUR    ,KSHAVE
C
         DIMENSION
     R   AFBS(15)  ,AFBS2(15) ,AP(15)   ,AP2(15)  ,AT(15)   ,
     R   AT2(15)   ,AT3(15)   ,AT4(15)

