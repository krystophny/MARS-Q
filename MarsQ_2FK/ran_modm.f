      MODULE RAN_MODM
      IMPLICIT NONE
      CONTAINS
C   ---------------
      FUNCTION RAN()
      IMPLICIT NONE
      INTEGER, SAVE :: FLAG = 0
      DOUBLE PRECISION :: RAN
      IF (FLAG == 0) THEN
         CALL RANDOM_SEED()
         FLAG = 1
      END IF
      CALL RANDOM_NUMBER(RAN)
      END FUNCTION RAN

C   ---------------
      FUNCTION RAN_NORMAL(MEAN, SIGMA)
      IMPLICIT NONE
      INTEGER :: FLAG
      DOUBLE PRECISION, PARAMETER :: PI = 3.141592653589793239
      DOUBLE PRECISION :: U1,U2,Y1,Y2,RAN_NORMAL,MEAN,SIGMA
      SAVE FLAG
      DATA FLAG /0/
      U1 = RAN()
      U2 = RAN()
      IF (FLAG. EQ. 0) THEN
         Y1 = SQRT(-2.0D0*LOG(U1))*COS(2.0D0*PI*U2)
         RAN_NORMAL = MEAN + SIGMA*Y1
         FLAG = 1
      ELSE
         Y2 = SQRT(-2.0D0*LOG(U1))*COS(2.0D0*PI*U2)
         RAN_NORMAL = MEAN + SIGMA*Y2
         FLAG = 0
      ENDIF
      END FUNCTION RAN_NORMAL
      END MODULE RAN_MODM

C   --------------------
C      PROGRAM TEST_RAN_MOD
C      USE RAN_MOD
C      IMPLICIT NONE
C      INTEGER, PARAMETER :: N = 10000
C      REAL(KIND = 8) :: A(N)
C      INTEGER :: I
C      OPEN (12, FLIE = '')
C      DO I = 1, N
C         A(I) = NORMAL(5.0D0, 2.0D0)
C         WRITE (12, *) A(I)
C      END DO
C      CLOSE (12)
C      END PROGRAM TEST_RAN_MOD 
