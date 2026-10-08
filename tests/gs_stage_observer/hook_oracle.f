      PROGRAM HOOK_ORACLE
      INCLUDE 'DECLAR.inc'
      INCLUDE 'COMDIM.inc'
      INCLUDE 'COMBLA.inc'
      INCLUDE 'COMBND.inc'
      INCLUDE 'COMCON.inc'
      INCLUDE 'COMNUM.inc'
      INCLUDE 'COMINT.inc'
      INCLUDE 'COMSOL.inc'
      INCLUDE 'COMPHY.inc'
      INCLUDE 'COMESH.inc'
      CHARACTER*1024 OUT_DIR
      DIMENSION SIG_READ(3),IDIMS(3)
      CALL GET_COMMAND_ARGUMENT(1,OUT_DIR)
      NS=2
      NT=4
      NS1=3
      NT1=5
      N4NSNT=48
      NBAND=3
      NBPS=3
      R0=1D0
      RZ0=0D0
      RMAG=R0
      RZMAG=RZ0
      SPSIM=-0.01D0
      R0EXP=6.2D0
      B0EXP=5.3D0
      RELAX=0.2D0
      A(1:NBAND,1:N4NSNT)=0D0
      CS(1:3)=(/0D0,0.99D0,0.999D0/)
      CSIG(1:3)=(/0D0,0.23D0,1D0/)
      DO I=1,NT1
        CT(I)=DFLOAT(I-1)
      ENDDO
      DO I=1,NBPS
        TETBPS(I,1)=DFLOAT(I-1)
        RRBPS(I,1)=1.1D0
        RZBPS(I,1)=0D0
        D2RBPS(I,1)=0D0
        D2ZBPS(I,1)=0D0
      ENDDO
      CALL GS_CAPTURE_MATRIX
      OPEN(NEWUNIT=IU,FILE=TRIM(OUT_DIR)//'/g0001_chart.dat',
     & FORM='UNFORMATTED')
      READ(IU) IDIMS
      READ(IU) SIG_READ
      CLOSE(IU)
      IF (ANY(IDIMS.NE.(/NS,NT,NBPS/))) STOP 1
      IF (ANY(SIG_READ.NE.(/0D0,0.23D0,1D0/))) THEN
        WRITE(*,*) 'WRONG RADIAL GRID: ',SIG_READ
        STOP 2
      ENDIF
      END
