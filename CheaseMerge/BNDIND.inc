C*COMDECK BNDIND
C ----------------------------------------------------------------------
C --     STATEMENT FUNCTION FOR BAND MATRIX ADRESSES                  --
C --                         11.07.88            HL        CRPP       --
C --                                                                  --
C --     INDCOL AND INDROW ARE COMPUTING BAND MATRIX ADRESSES STARTING--
C --     FROM SQUARE MATRIX ADRESSES. IF THE ADRESS IS LOWER-DIAGONAL,--
C --     THESE FUNCTIONS FIND OUT AUTOMATICALY THE SYMETRIC UPPER-    --
C --     DIAGONAL ADRESS.                                             --
C --                                                                  --
C ----------------------------------------------------------------------
C
         INDCOL(ICOL,IROW) = SIGN(1,ICOL - IROW) * (ICOL - IROW) + 1
         INDROW(ICOL,IROW) = ((IROW + ICOL) + (IROW - ICOL) *
     *                        SIGN(1,ICOL - IROW)) / 2
C
