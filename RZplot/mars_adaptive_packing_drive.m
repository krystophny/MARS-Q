%radial grid number in new adaptive packed equilibrium
new_NRP1 = 301
%path and filename of equilibrium used in the previous MARS run
old_equlocate='OUTRMAR'
%location of eigenfunction in the previous MARS run
eigenfile_dir=''
%path and file of high resolution equilibrium
highresol_equlocate='OUTRMAR_HIGHRESOL'
%path and file of new requilibruim to output
new_equlocate='OUTRMAR'

mars_adaptive_packing_routines (new_NRP1,old_equlocate,highresol_equlocate,new_equlocate, eigenfile_dir);