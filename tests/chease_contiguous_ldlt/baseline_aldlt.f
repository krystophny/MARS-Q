         SUBROUTINE BASELINE_ALDLT(A,EPS,N,M,MP,NSING)
C        ------------------------------------
C
C     DECOMPOSE A=L*D*LT                                             
C                                                                    
C     VERSION 1C           13.9.74     RALF GRUBER    CRPP LAUSANNE  
C                                                                    
C     A IS A BAND MATRIX WITH HALF WIDTH M AND LENGTH N              
C     L CONTAINS 1 IN THE DIAGONAL                                   
C     AS OUTPUT D REPLACES THE DIAGONAL OF A AND                     
C     LT WITHOUT ITS DIAGONAL THE REST OF A                          
C     ALL CALCULATIONS ARE PERFORMED IN A                            
C     NSING = -1 WHEN A IS SINGULAR                                  
C
C
         INCLUDE 'DECLAR.inc'
         DIMENSION
     R   A(N*MP)
C
C     INITIALIZE
C
         M1  = MP - 1
         IKD = 0
         AD  = ABS(A(1)) * EPS
C
C     SCAN OVER THE WHOLE LENGTH OF A
C
         DO  4  JIB=2,N
            DIAG=A(IKD+1)
C
C     TEST FOR ZERO PIVOT
C
            IF (ABS(DIAG) .LT. AD) THEN
               NSING = -1
               RETURN
            ENDIF
C
C     RESTRICTION OF LOOP FOR NOT EXCEEDING BAND MATRIX
C
            LOPBND = M
            I1     = N - JIB + 2
C
            IF (I1 .LT. M) LOPBND = I1
C
C     DIAGONAL ELEMENT BEFORE GAUSS ELIMINATION
C
            IJ = IKD + MP + 1
            AD = ABS(A(IJ)) * EPS
C
C     SETS THE ROW OF THE TRANSPOSED LEFT HAND SIDE MATRIX LT
C
            DO 3 JJB=2,LOPBND
               ITOP    = IKD + JJB
               TOP     = A(ITOP)
               A(ITOP) = A(ITOP) / DIAG
C
C     GAUSS RECTANGULAR RULE GOING DOWNWARDS
C
               CALL SAXPY(JJB-1,-TOP,A(IKD+2),1,A(ITOP+M1),M1)
   3        CONTINUE
         IKD = IKD + MP
   4     CONTINUE
C
C     LAST DIAGONAL ELEMENT
C
         IKD = (N - 1) * MP + 1
         IJ  = IKD - MP
C
         IF (ABS(A(IKD)) .LT. ABS(A(IJ))*EPS) THEN
            NSING = -1
         ENDIF
C
         RETURN
         END
