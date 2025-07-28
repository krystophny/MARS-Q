C*COMDECK QUAQDQ
C ----------------------------------------------------------------------
C --     STATEMENT FUNCTION FOR QUADRATIC INTERPOLATION               --
C --                         19.01.87            AR        CRPP       --
C --                                                                  --
C -- QUADRATIC INTERPOLATION OF A FUNCTION F(X)                       --
C -- THE FIVE PARAMETERS X1,F1,P1,X2,F2    ARE DEFINED AS FOLLOWS:    --
C -- F(X1) = F1 , DF/DX(X1) = P1 , F(X2) = F2                         --
C ----------------------------------------------------------------------
C
         FD2(X1,F1,P1,X2,F2) = ((F2-F1)/(X2-X1) - P1) / (X2-X1)
         FD1(X1,F1,P1,X2,F2) = P1 - 2.*X1*FD2(X1,F1,P1,X2,F2)
         FD0(X1,F1,P1,X2,F2) = F1 - X1*(X1*FD2(X1,F1,P1,X2,F2) +
     +                                     FD1(X1,F1,P1,X2,F2))
C ----------------------------------------------------------------------
C -- FQDQ0 GIVES THE VALUE OF THE FUNCTION AT POINT PX                --
C -- FQDQ0(......,PX) = F(PX)                                         --
C ----------------------------------------------------------------------
         FQDQ0(X1,F1,P1,X2,F2,PX) = FD0(X1,F1,P1,X2,F2) +
     F                              PX * (FD1(X1,F1,P1,X2,F2) +
     F                                    PX * FD2(X1,F1,P1,X2,F2))
C ----------------------------------------------------------------------
C -- FQDQ1 GIVES THE VALUE OF THE FIRST DERIVATIVE OF F(X) AT PX      --
C -- FQDQ1(......,PX) = DF/DX (PX)                                    --
C ----------------------------------------------------------------------
         FQDQ1(X1,F1,P1,X2,F2,PX) = FD1(X1,F1,P1,X2,F2) +
     F                              2.* PX * FD2(X1,F1,P1,X2,F2)
C ----------------------------------------------------------------------
C -- FQDQ2 GIVES THE VALUE OF THE SECOND DERIVATIVE OF F(X) AT PX     --
C -- FQDQ2(......,PX) = D2F/DX2 (PX)                                  --
C ----------------------------------------------------------------------
         FQDQ2(X1,F1,P1,X2,F2) = 2. * FD2(X1,F1,P1,X2,F2)
