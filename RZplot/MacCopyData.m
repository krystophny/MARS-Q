mach='AUG';
numb='33133';
time='3000';
casa='';
%cass=['a=' casa '_RMP_EF'];
cass=['PLS_WF_N'];

shot=[mach numb];
if mach=='AUG'
   n=2; mk=-29:29;
elseif mach=='D3D'
   n=3; mk=-33:33;
end 
   


eval(['copyfile ../' shot '/Data/RMZM_F_EQAC_' shot '.' time  ' ../Temp/RMZM_F_GEOM']);
eval(['copyfile ../' shot '/Data/RMZM_F_PEST_' shot '.' time  ' ../Temp/RMZM_F_PEST']);
eval(['copyfile ../' shot '/Data/PROFEQ_' shot '.' time  ' ../Temp/PROFEQ.OUT']);
eval(['copyfile ../' shot '/Data_New/PROFROT_' shot '.' time 'N.IN ../Temp/PROFROT.IN']);
eval(['copyfile ../' shot '/Data_New/PROFWE_' shot '.' time 'N.IN ../Temp/PROFWE.IN']);
eval(['copyfile ../' shot '/Data_New/BPLASMA_EQAC_' shot '.' time  '_' cass ' ../Temp/BPLASMA.OUT']);
if 1==1
   eval(['copyfile ../' shot '/Data_New/XPLASMA_EQAC_' shot '.' time  '_' cass ' ../Temp/XPLASMA.OUT']);
end

if 1==1
MacRfaCtBn2
save b1res.txt BnPEST_RS -ascii 
eval(['!mv b1res.txt ../' shot '/Data_New/B1res_' shot '.' time '_' cass])
eval(['!mv profq.fig ../' shot '/Data_New/ProfQ_' shot '.' time  '_' cass '.fig'])
eval(['!mv profw.fig ../' shot '/Data_New/ProfRot_' shot '.' time '_' cass '.fig'])
eval(['!mv b1m2d.fig ../' shot '/Data_New/B1m2D_' shot '.' time '_' cass '.fig'])
eval(['!mv b1m1d.fig ../' shot '/Data_New/B1m1D_' shot '.' time '_' cass '.fig'])
end

if 1==1
   MacRfaCtVn2
  eval(['!mv XnSurfGeom.txt ../' shot '/Data_New/XnSurfGeom_' shot '.' time '_' cass])
  eval(['!mv x1m1d.fig ../' shot '/Data_New/X1m1D_' shot '.' time '_' cass '.fig'])
  eval(['!mv xn2d.fig ../' shot '/Data_New/Xn2D_' shot '.' time '_' cass '.fig'])
end
