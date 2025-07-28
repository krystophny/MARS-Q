C*COMDECK QUAQQQ
C ----------------------------------------------------------------------
C --     STATEMENT FUNCTION FOR QUADRATIC INTERPOLATION               --
C --                         19.01.87            AR        CRPP       --
C --                                                                  --
C -- QUADRATIC INTERPOLATION OF A FUNCTION F(X)                       --
C -- THE SIX PARAMETERS A1,A2,A3,B1,B2,B3 ARE DEFINED AS FOLLOWS:     --
C -- F(B1) = A1 , F(B2) = A2 , F(B3) = A3                             --
C ----------------------------------------------------------------------
C
         FB2(A1,A2,A3,B1,B2,B3) =
     F               ((A1-A2)/(B1-B2)-(A1-A3)/(B1-B3))/(B2-B3)
         FB1(A1,A2,A3,B1,B2,B3) = ((A1-A2)/(B1-B2))-
     F         FB2(A1,A2,A3,B1,B2,B3)*(B1+B2)
         FB0(A1,A2,A3,B1,B2,B3) = A1-FB1(A1,A2,A3,B1,B2,B3)*B1
     F         -FB2(A1,A2,A3,B1,B2,B3)*B1*B1
C ----------------------------------------------------------------------
C -- FQQQ0 GIVES THE VALUE OF THE FUNCTION AT THE POINT PX            --
C -- FQQQ0(......,PX) = F(PX)                                         --
C ----------------------------------------------------------------------
         FQQQ0(A1,A2,A3,B1,B2,B3,PX) = FB0(A1,A2,A3,B1,B2,B3) +
     F                                 PX * (FB1(A1,A2,A3,B1,B2,B3) +
     F                                 PX * FB2(A1,A2,A3,B1,B2,B3))
C ----------------------------------------------------------------------
C -- FQQQ1 GIVES THE VALUE OF THE FIRST DERIVATIVE OF F(X) AT PX      --
C -- FQQQ1(......,PX) = DF/DX (PX)                                    --
C ----------------------------------------------------------------------
         FQQQ1(A1,A2,A3,B1,B2,B3,PX) = FB1(A1,A2,A3,B1,B2,B3) +
     F                                 2. * PX * FB2(A1,A2,A3,B1,B2,B3)
C ----------------------------------------------------------------------
C -- FQQQ2 GIVES THE VALUE OF THE SECOND DERIVATIVE OF F(X) AT PX     --
C -- FQQQ2(......,PX) = D2F/DX2 (PX)                                  --
C ----------------------------------------------------------------------
         FQQQ2(A1,A2,A3,B1,B2,B3) = 2. * FB2(A1,A2,A3,B1,B2,B3)
