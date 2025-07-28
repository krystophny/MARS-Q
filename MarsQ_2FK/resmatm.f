      MODULE RESMATM
C
C..........  ORGANISATION OF MATRICES AND VECTORS IN MARS ..............
C......   (MATRICES A & H NUMBERED DIFFERENTLY THAN IN PAMS) ...........
C
C        B1 E1 C1                      X1
C        F1 D1 G1                      Y1
C        A2 H2 B2 E2 C2                X2
C              F2 D2 G2                Y2
C              A3 H3 B3 E3 C3          X3
C
c-----------------------------------------------------------------------
C
      COMPLEX*16,DIMENSION(:,:,:),ALLOCATABLE::ASUBM,BSUBM,CSUBM,
     &     DSUBM,ESUBM,FSUBM,GSUBM,HSUBM
C      COMPLEX*16          ASUBM(MXMAT*MXMAT*NCRAY),
C     &                 BSUBM(MXMAT*MXMAT*NCRAY),
C     &                 CSUBM(MXMAT*MXMAT*NCRAY),
C     &                 DSUBM(MYMAT*MYMAT*NCRAY),
C     &                 ESUBM(MXMAT*MYMAT*NCRAY),
C     &                 FSUBM(MYMAT*MXMAT*NCRAY),
C     &                 GSUBM(MYMAT*MXMAT*NCRAY),
C     &                 HSUBM(MXMAT*MYMAT*NCRAY)
CC
C      COMMON /MATRIX/ ASUBM,BSUBM,CSUBM,DSUBM,ESUBM,FSUBM,GSUBM,HSUBM
C
C-----------------------------------------------------------------------
C
      COMPLEX*16,DIMENSION(:,:),ALLOCATABLE::DX,DY,X,Y
C      COMPLEX*16   DX(MXMAT,NCRAY),DY(MYMAT,NCRAY)
C      COMPLEX*16    X(MXMAT,NCRAY), Y(MYMAT,NCRAY)
C
C-----------------------------------------------------------------------
C
      INTEGER     MD,MDY
C      PARAMETER   (MD=MXMAT,MDY=MYMAT)
C      INTEGER IWORK
C      PARAMETER (IWORK = 2*( 3*NCRAY*(MD+MDY) + 10*MD*MD+2*MD*MDY
C     &                           + 4*MDY*MDY + 10*MD + 7*MDY ))
C
         REAL*8,DIMENSION(:),ALLOCATABLE::WORK
C         COMMON/COWORK/ WORK,DX,DY,X,Y
C-----------------------------------------------------------------------
         END MODULE RESMATM




