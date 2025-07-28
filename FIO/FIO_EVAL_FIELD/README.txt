1.Compile fio*.f using: make ==> fio_test.x
  Need to change source code fio_test.f to specify (R,Z) location

2. run MARS-F (new version) 
   change input variables from OUTOPT: ORMIN=0.8,ORMAX=2.6,OZMIN=-1.6,OZMAX=1.6,NORR=91,NOZZ=161  ==>
   SCHIMESH_RECTRZ0001.OUT ==> copy to SCHIMESH_RECTRZ.IN
   BPLASMA_MARSF0001.OUT   ==> copy to BPLASMA_MARSF_1.IN

3. run fio_test.x ==> record deltaB value
 
