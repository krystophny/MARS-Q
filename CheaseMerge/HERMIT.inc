C*COMDECK HERMIT
C ----------------------------------------------------------------------
C --     STATEMENT FUNCTION FOR CUBIC HERMIT BASIS FUNCTIONS ,        --
C --     FIRST AND SECOND DERIVATIVES OF CUBIC HERMIT BASIS FUBCTIONS --
C --                                                                  --
C --                         11.07.88            HL        CRPP       --
C --                                                                  --
C ----------------------------------------------------------------------
C --                                                                  --
C -- 1.  ONE DIMENSION HERMIT CUBIC BASIS FUNCTIONS                   --
C --                                                                  --
C ----------------------------------------------------------------------
C
         FN1(X,XJ,XJ1) = (X - XJ1) * (X - XJ1) * (XJ1 - 3*XJ + 2*X) /
     /                   ((XJ1 - XJ) * (XJ1 - XJ) * (XJ1 - XJ))
         FN3(X,XJ,XJ1) = (X - XJ) * (X - XJ) * (3*XJ1 - XJ - 2*X) /
     /                   ((XJ1 - XJ) * (XJ1 - XJ) * (XJ1 - XJ))
         FN2(X,XJ,XJ1) = (X - XJ1) * (X - XJ1) * (X - XJ) /
     /                   ((XJ1 - XJ) * (XJ1 - XJ))
         FN4(X,XJ,XJ1) = (X - XJ) * (X - XJ) * (X - XJ1) /
     /                   ((XJ1 - XJ) * (XJ1 - XJ))
C
C ----------------------------------------------------------------------
C --                                                                  --
C -- 2. FIRST DERIVATIVE OF ONE DIMENSION CUBIC HERMIT BASIS FUNCTIONS--
C --                                                                  --
C ----------------------------------------------------------------------
C
         DFN1(X,XJ,XJ1) = 6 * (X - XJ1) * (X - XJ) /
     /                    ((XJ1 - XJ) * (XJ1 - XJ) * (XJ1 - XJ))
         DFN3(X,XJ,XJ1) = 6 * (XJ1 - X) * (X - XJ) /
     /                    ((XJ1 - XJ) * (XJ1 - XJ) * (XJ1 - XJ))
         DFN2(X,XJ,XJ1) = (X - XJ1) * (3 * X - 2 * XJ - XJ1) /
     /                    ((XJ1 - XJ) * (XJ1 - XJ))
         DFN4(X,XJ,XJ1) = (X - XJ) * (3 * X - XJ - 2 * XJ1) /
     /                    ((XJ1 - XJ) * (XJ1 - XJ))
C
C ----------------------------------------------------------------------
C --                                                                  --
C -- 3. SECOND DERIVATIVE OF ONE DIMENSION CUBIC HERMIT BASIS         --
C --    FUNCTIONS                                                     --
C --                                                                  --
C ----------------------------------------------------------------------
C
         D2FN1(X,XJ,XJ1) = 6 * (2 * X - XJ - XJ1) /
     /                     ((XJ1 - XJ) * (XJ1 - XJ) * (XJ1 - XJ))
         D2FN3(X,XJ,XJ1) = 6 * (XJ + XJ1 - 2 * X) /
     /                     ((XJ1 - XJ) * (XJ1 - XJ) * (XJ1 - XJ))
         D2FN2(X,XJ,XJ1) = 2 * (3 * X - XJ - 2 * XJ1) /
     /                     ((XJ1 - XJ) * (XJ1 - XJ))
         D2FN4(X,XJ,XJ1) = 2 * (3 * X - 2 * XJ - XJ1) /
     /                     ((XJ1 - XJ) * (XJ1 - XJ))
C
