C*COMDECK CUCDCD
C ----------------------------------------------------------------------
C --     STATEMENT FUNCTION FOR CUBIC INTERPOLATION                   --
C --                         19.01.87            AR        CRPP       --
C --                                                                  --
C -- CUBIC INTERPOLATION OF A FUNCTION F(X)                           --
C -- THE SIX ARGUMENTS X1,F1,P1,X2,F2,P2 ARE DEFINED AS FOLLOWS:      --
C -- F(X1) = F1 , F(X2) = F2 , DF/DX(X1) = P1 , DF/DX(X2) = P2        --
C ----------------------------------------------------------------------
C
         FC3(X1,F1,P1,X2,F2,P2) =
     =      (2. * (F2 - F1) / (X1 - X2) + (P1 + P2)) / 
     /      ((X1 - X2) * (X1 - X2))
         FC2(X1,F1,P1,X2,F2,P2) =
     =      (3. * (X1 + X2) * (F1 - F2) / (X1 - X2) - 
     *       P1 * (X1 + 2. * X2) - P2 * (X2 + 2. * X1)) /
     /      ((X1 - X2) * (X1 - X2))
         FC1(X1,F1,P1,X2,F2,P2) =
     =      (6. * X1 * X2 * (F2 - F1) / (X1 - X2) + 
     *       X2 * P1 * (2 * X1 + X2) + X1 * P2 * (X1 + 2. * X2)) /
     /      ((X1 - X2) * (X1 - X2))
         FC0(X1,F1,P1,X2,F2,P2) =
     =      (F1 * X2**2 + F2 * X1**2 - X1 * X2 * (X2 * P1 + X1 * P2) +
     +       2. * X1 * X2 * (F1 * X2 - F2 * X1) / (X1 - X2)) /
     /      ((X1 - X2) * (X1 - X2))
C ----------------------------------------------------------------------
C -- FCDCD0 GIVES THE VALUE OF THE FUNCTION AT POINT PX               --
C -- FCDCD0(......,PX) = F(PX)                                        --
C ----------------------------------------------------------------------
         FCDCD0(X1,F1,P1,X2,F2,P2,PX) =
     =              FC0(X1,F1,P1,X2,F2,P2) +
     +              PX * (FC1(X1,F1,P1,X2,F2,P2) +
     +                    PX * (FC2(X1,F1,P1,X2,F2,P2) +
     +                          PX * FC3(X1,F1,P1,X2,F2,P2)))
C ----------------------------------------------------------------------
C -- FCDCD1 GIVES THE VALUE OF THE FIRST DERIVATIVE OF F(X) AT PX:    --
C -- FCDCD1(......,PX) = DF/DX (PX)                                   --
C ----------------------------------------------------------------------
         FCDCD1(X1,F1,P1,X2,F2,P2,PX) =
     =              FC1(X1,F1,P1,X2,F2,P2) +
     +              PX * (2. * FC2(X1,F1,P1,X2,F2,P2) +
     +                    3. * PX * FC3(X1,F1,P1,X2,F2,P2))
C ----------------------------------------------------------------------
C -- FCDCD2 GIVES THE VALUE OF THE SECOND DERIVATIVE OF F(X) AT PX:   --
C -- FCDCD2(......,PX) = D2F/DX2 (PX)                                 --
C ----------------------------------------------------------------------
         FCDCD2(X1,F1,P1,X2,F2,P2,PX) =
     =             2. * FC2(X1,F1,P1,X2,F2,P2) +
     +             6. * FC3(X1,F1,P1,X2,F2,P2) * PX
C ----------------------------------------------------------------------
C -- FCDCD3 GIVES THE VALUE OF THE THIRD DERIVATIVE OF F(X) AT PX:    --
C -- FCDCD3(......,PX) = D3F/DX3 (PX)                                 --
C ----------------------------------------------------------------------
         FCDCD3(X1,F1,P1,X2,F2,P2,PX) =
     =                      6. * FC3(X1,F1,P1,X2,F2,P2)
