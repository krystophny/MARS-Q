C*COMDECK CUCCCC
C ----------------------------------------------------------------------
C --     STATEMENT FUNCTION FOR CUBIC INTERPOLATION                   --
C --                         23.04.88            AR        CRPP       --
C --                                                                  --
C -- CUBIC INTERPOLATION OF A FUNCTION F(X)                           --
C -- THE EIGHT ARGUMENTS A1,A2,A3,A4,B1,B2,B3,B4 ARE DEFINED BY:      --
C -- F(B1) = A1 , F(B2) = A2 , F(B3) = A3 , F(B4) = A4                --
C ----------------------------------------------------------------------
C
         FA3(A1,A2,A3,A4,B1,B2,B3,B4) =
     F        (A1-A2) / ((B1-B2)*(B2-B4)*(B2-B3)) +
     F        (A1-A3) / ((B4-B3)*(B3-B1)*(B3-B2)) +
     F        (A1-A4) / ((B1-B4)*(B2-B4)*(B3-B4))
         FA2(A1,A2,A3,A4,B1,B2,B3,B4) =
     F        (A1-A2) / ((B2-B1)*(B3-B2)) +
     F        (A3-A1) / ((B3-B1)*(B3-B2)) -
     F        (B1+B2+B3) * FA3(A1,A2,A3,A4,B1,B2,B3,B4)
         FA1(A1,A2,A3,A4,B1,B2,B3,B4) =
     F        (A1-A2) / (B1-B2) -
     F        (B1+B2) * FA2(A1,A2,A3,A4,B1,B2,B3,B4) -
     F        (B1*B1+B1*B2+B2*B2) * FA3(A1,A2,A3,A4,B1,B2,B3,B4)
         FA0(A1,A2,A3,A4,B1,B2,B3,B4) =
     F        A1 -
     F        B1 * (FA1(A1,A2,A3,A4,B1,B2,B3,B4) +
     F              B1 * (FA2(A1,A2,A3,A4,B1,B2,B3,B4) +
     F                    B1 * FA3(A1,A2,A3,A4,B1,B2,B3,B4)))
C ----------------------------------------------------------------------
C -- FCCCC0 GIVES THE VALUE OF THE FUNCTION AT POINT PX:              --
C -- FCCCC0(......,PX) = F(PX)                                        --
C ----------------------------------------------------------------------
        FCCCC0(A1,A2,A3,A4,B1,B2,B3,B4,PX) =
     F              FA0(A1,A2,A3,A4,B1,B2,B3,B4) +
     F              PX * (FA1(A1,A2,A3,A4,B1,B2,B3,B4) +
     F                    PX * (FA2(A1,A2,A3,A4,B1,B2,B3,B4) +
     F                          PX * FA3(A1,A2,A3,A4,B1,B2,B3,B4)))
C ----------------------------------------------------------------------
C -- FCCCC1 GIVES THE VALUE OF THE FIRST DERIVATIVE OF F(X) AT PX:    --
C -- FCCCC1(......,PX) = DF/DX (PX)                                   --
C ----------------------------------------------------------------------
        FCCCC1(A1,A2,A3,A4,B1,B2,B3,B4,PX) =
     F              FA1(A1,A2,A3,A4,B1,B2,B3,B4) +
     F              PX * (2. * FA2(A1,A2,A3,A4,B1,B2,B3,B4) +
     F                    3. * PX * FA3(A1,A2,A3,A4,B1,B2,B3,B4))
C ----------------------------------------------------------------------
C -- FCCCC2 GIVES THE VALUE OF THE SECOND DERIVATIVE OF F(X) AT PX:   --
C -- FCCCC2(......,PX) = D2F/DX2 (PX)                                 --
C ----------------------------------------------------------------------
         FCCCC2(A1,A2,A3,A4,B1,B2,B3,B4,PX) =
     F             2. * FA2(A1,A2,A3,A4,B1,B2,B3,B4) +
     F             6. * FA3(A1,A2,A3,A4,B1,B2,B3,B4) * PX
C ----------------------------------------------------------------------
C -- FCCCC3 GIVES THE VALUE OF THE THIRD DERIVATIVE OF F(X) AT PX:     -
C -- FCCCC3(......,PX) = D3F/DX3 (PX)                                  -
C ----------------------------------------------------------------------
         FCCCC3(A1,A2,A3,A4,B1,B2,B3,B4,PX) =
     F                      6. * FA3(A1,A2,A3,A4,B1,B2,B3,B4)
