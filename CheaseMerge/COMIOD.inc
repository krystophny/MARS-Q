C*COMDECK COMIOD
C
***********************************************************************
*                                                                     *
*     C.09         CONTROL I/O                                        *
*                                                                     *
*     DISK CHANNELS :                                                 *
*     ---------------                                                 *
*     FILE CONTENT                       UNIT #                       *
*                                                                     *
*     EQ'S FOR NOVA-W AND PEST         : INP1
*     STORE EQ'S FOR ERATO             : MEQ                          *
*     STORE PLOT QUANTITIES FOR ERATO  : NDES                         *
*     STORE VACUUM QUANTITIES FOR MARS : NETVAC                       *
*     READ OLD EQUILIBRIUM             : NIN                          *
*     STORE MARS EQ'S                  : NO                           *
*     STORE QUANTITIES FOR MARS-FLOW   : NOI                          *
*     STORE QUANTITIES FOR MARS-FLOW   : NO3                          *
*     STORE COMPUTED EQUILIBRIUM       : NOUT                         *
*     EQ'S FOR PENN                    : NOUT                         *
*     STORE OUTPUT                     : NPRNT                        *
*     STORE BOUNDARY OF SHAVED                                        *
*     EQUILIBRIUM                      : NRZPEL                       *
*     STORE NAMELIST AND LABELS        : NSAVE                        *
*     STORE QUANTITIES FOR EQUILIBRIUM                                *
*           PLOTS                      : NUPLO                        *
*     STORE VACUUM QUANTITIES FOR ERATO: NVAC                         *
*     READ EXPERIMENTAL EQUILIBRIUM    : NXIN                         *
*     RESTORE EXPERIMENTAL EQUILIBRIUM : NXOUT                        *
*     STORE EQ'S FOR XTOR              : NXTOR                        *
*                                                                     *
***********************************************************************
C
         COMMON /COMIOD/
     I       MEQ    ,NDES   ,NETVAC ,NIN    ,NO     ,INP1,
     I       NOI    ,NO3    ,NOUT   ,NPENN  ,NPRNT  ,NRZPEL ,
     I       NSAVE  ,NUPLO  ,NVAC   ,NXIN   ,NXOUT  ,NXTOR  ,
     I       NRMAR  ,NO4    ,NO5    ,NRDCON

