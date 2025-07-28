%Combine vacuum field from Aalto and MARS-F computed field
%BPLASMA_RECTRZ_BIOT, to obtain the total field which is
%valid in the whole rectangular (R,Z) domain

%close all, 
clear all,
kplot2d = 0;
kdo = 2;  %1: compare vacuum field between Aalto and MARS-F computations
              %2: combine Aalto vacuum field with MARS-F BIOT-Savart field     
              %3: read MARS-F compute response field by ELM coils     
%SDIR  = 'C:\Users\Yueqiang\Documents\LIU\Temp\Temp\F4E_Contract_15001\Work\';
SDIR  = '/home/yliu/Temp/';
SSAV = '/home/yliu/F4E_2014/DataRes/10100/n=1/9MAFTFITBMELM/';

if kdo==1
eval(['!cp ' SDIR 'BPLASMA_RECTRZ_VAC ' SDIR 'BPLASMA_RECTRZ']);
LCR = 'b';
SLW = 2;
KASCOT = 0;
MacReadRectRZ

LCR = 'r';
SLW = 2;
KASCOT = 1;
MacReadRectRZ
end

if kdo==2
%eval(['!cp ' SDIR 'BPLASMA_RECTRZ_BIOT ' SDIR 'BPLASMA_RECTRZ']);
eval(['!cp ' SDIR 'BPLASMA_RECTRZ_BIOT_KIN ' SDIR 'BPLASMA_RECTRZ']);
LCR = 'r';
SLW = 2;
KASCOT = 0;
MacReadRectRZ
BnBIOT = Bn;
   
LCR = 'b';
SLW = 2;
KASCOT = 1;
MacReadRectRZ
BnVAC = Bn;

Bn.Br = BnVAC.Br + BnBIOT.Br;
Bn.Bz = BnVAC.Bz + BnBIOT.Bz;
Bn.Bphi = BnVAC.Bphi + BnBIOT.Bphi;
save BnTOT Bn
eval(['!mv BnTOT.mat ' SSAV '.']);
LCR = 'g';
SLW = 4;
KASCOT = 2;
MacReadRectRZ

%eval(['!cp ' SDIR 'BPLASMA_RECTRZ_PLS ' SDIR 'BPLASMA_RECTRZ']);
eval(['!cp ' SDIR 'BPLASMA_RECTRZ_PLS_KIN ' SDIR 'BPLASMA_RECTRZ']);
LCR = 'k';
SLW = 4;
KASCOT = 0;
MacReadRectRZ

figure(14)
print -depsc BRvsRCmp
%eval(['!mv BRvsRCmp.eps ' SSAV '.'])
end

if kdo==3
eval(['!cp ' SDIR 'BPLASMA_RECTRZ_ELM ' SDIR 'BPLASMA_RECTRZ']);
LCR = 'r';
SLW = 4;
KASCOT = 0;
MacReadRectRZ
save BnELM Bn
eval(['!mv BnELM.mat ' SSAV '.']);
end
