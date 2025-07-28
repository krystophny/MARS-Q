1. Compile fio*.f using: make ==> fio_main.x
   Change source code fio_main.f to specify new (R,PHI,Z) location

2. fio_main.x reads in the following input files as direct output from MARS-F
   SCHIMESH_RECTRZ.IN
   BPLASMA_MARSF_1.IN

3. run fio_main.x ==> output (B_R,B_PHI,B_Z) AND ALL FIRST ORDER DERIVATIVES W.R.T. (R,PHI,Z) 
 
