global SDIR 

kcase = 1;  %1: MAST, 2: D3D139571, 3: ANATOR

if kcase == 1

SDIR = '/.automount/funsrv1/root/home/yliu/MAST_ALL/';  

shot = '24460'; flow='_IFLOW'; coil='_ODD'; q95n='_3'; wf=''; kclr = 1; Nclr=2;

eval(['!cp -f /scratch/yliu/Mast/MAST_ALL/RMZM_F_PEST_' shot q95n ' ' SDIR 'RMZM_F_PEST']);
eval(['!cp -f /scratch/yliu/Mast/MAST_ALL/RMZM_F_EQAC_' shot q95n ' ' SDIR 'RMZM_F_EQAC']);
eval(['!cp -f ' SDIR 'PROFEQ_MAST' shot q95n ' ' SDIR 'PROFEQ_PEST']);
%eval(['!cp -f ' SDIR 'BnMat_VAC_MAST' shot coil '.mat ' SDIR 'BnMat_VAC.mat']);
%eval(['!cp -f ' SDIR 'BnMat' flow coil '_MAST' shot '.mat ' SDIR 'BnMat.mat']);

%eval(['!cp -f ' SDIR 'RMZM_F_EQAC ' SDIR 'RMZM_F']);
eval(['!cp -f /scratch/yliu/Mast/MAST_ALL/VPLASMA' flow coil '_' shot q95n wf ' ' SDIR 'VPLASMA']);
%eval(['!cp -f ' SDIR 'XnSurf_' shot '_5' flow coil '.txt  ' SDIR 'XnSurf.txt']);

elseif kcase == 2

SDIR = '/.automount/funsrv1/root/home/yliu/D3D139571/Rfa/';  

%S1='_CFBAL1.0'; S2='_ROTE4.5'; S3='_ETA2e-7'; kclr = 2; Nclr=8;
S1='_JP3.0'; S2='_ROTE0.5'; S3=''; kclr = 8; Nclr=8;


eval(['!cp -f ' SDIR 'RMZM_F_PEST' S1 ' ' SDIR 'RMZM_F_PEST']);
eval(['!cp -f ' SDIR 'RMZM_F_EQAC' S1 ' ' SDIR 'RMZM_F_EQAC']);
eval(['!cp -f ' SDIR 'PROFEQ_PEST' S1 ' ' SDIR 'PROFEQ_PEST']);
%eval(['!cp -f ' SDIR 'BnMat_VAC' S1 '.mat ' SDIR 'BnMat_VAC.mat']);
eval(['!cp -f ' SDIR 'BnMat_RES' S1 S2 S3 '.mat ' SDIR 'BnMat.mat']);


elseif kcase == 3

SDIR = '/.automount/funsrv1/root/home/yliu/ANATOR/';  

%S1='_CFBAL1.0'; S2='_ROTE4.5'; S3='_ETA2e-7'; kclr = 2; Nclr=8;
S1='_EQ1'; S2=''; S3=''; kclr = 8; Nclr=8;


eval(['!cp -f ' SDIR 'RMZM_F' S1 '_PEST' ' ' SDIR 'RMZM_F_PEST']);
eval(['!cp -f ' SDIR 'RMZM_F' S1 '_EQAC' ' ' SDIR 'RMZM_F_EQAC']);
eval(['!cp -f ' SDIR 'PROFEQ' S1 '_PEST' ' ' SDIR 'PROFEQ_PEST']);
eval(['!cp -f ' SDIR 'BnMat_Ioddv.mat ' SDIR 'BnMat_VAC.mat']);
eval(['!cp -f ' SDIR 'BnMat_Iodd.mat ' SDIR 'BnMat.mat']);


end

