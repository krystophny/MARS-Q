% combine YPLASMA.OUT from U and L rows of RMP coils
% with varying coil phasing, where Y=SY='B','X',or 'P'

function MacScanCoilPhase

global RunB RunV RunP SDIR

SDIR = '/home/liuy/Work/MAST_RWM/Data/27112/';
SY   = 'B';
coil_phase = 120; %in [degree]

RunB = 0;
RunV = 0;
RunP = 0;
if strcmp(SY,'B'), RunB=1; end
if strcmp(SY,'X'), RunV=1; end
if strcmp(SY,'P'), RunP=1; end

%read in data from U- and L-rows
YU = load([SDIR 'RFA_NEW/' SY 'PLASMA_U_KIN_INCDPHI1']);
YL = load([SDIR 'RFA_NEW/' SY 'PLASMA_L_KIN_INCDPHI1']);
%YU = load([SDIR 'RFA_NEW/' SY 'PLASMA_U_FLD']);
%YL = load([SDIR 'RFA_NEW/' SY 'PLASMA_L_FLD']);

a = exp(i*coil_phase/180*pi);

Y = MacCombYPLASMA(YU,YL,a,SY);

save([SDIR SY 'PLASMA.OUT'],'Y','-ascii');
MacMainMASTU_B


function [YC]=MacCombYPLASMA(YU,YL,a,SY)

MSMAX = floor(YU(1,1));
NRP1  = floor(YU(1,2));

if strcmp(SY,'B'), NSTR=MSMAX+2; end
if strcmp(SY,'X'), NSTR=MSMAX+NRP1+2; end
if strcmp(SY,'P'), NSTR=MSMAX+2; end

Y1U = YU(NSTR:end,1) + YU(NSTR:end,2)*i;
Y2U = YU(NSTR:end,3) + YU(NSTR:end,4)*i;
Y3U = YU(NSTR:end,5) + YU(NSTR:end,6)*i;

Y1L = YL(NSTR:end,1) + YL(NSTR:end,2)*i;
Y2L = YL(NSTR:end,3) + YL(NSTR:end,4)*i;
Y3L = YL(NSTR:end,5) + YL(NSTR:end,6)*i;

Y1C = Y1U+Y1L*a;
Y2C = Y2U+Y2L*a;
Y3C = Y3U+Y3L*a;

if strcmp(SY,'B') | strcmp(SY,'X')
   YC = [YU(1:NSTR-1,:); real(Y1C) imag(Y1C) real(Y2C) imag(Y2C) real(Y3C) imag(Y3C)];
elseif strcmp(SY,'P')
   Y4U = YU(NSTR:end,7) + YU(NSTR:end,8)*i;
   Y5U = YU(NSTR:end,9) + YU(NSTR:end,10)*i;

   Y4L = YL(NSTR:end,7) + YL(NSTR:end,8)*i;
   Y5L = YL(NSTR:end,9) + YL(NSTR:end,10)*i;

   Y4C = Y4U+Y4L*a;
   Y5C = Y5U+Y5L*a;

   YC = [YU(1:NSTR-1,:); real(Y1C) imag(Y1C) real(Y2C) imag(Y2C) real(Y3C) imag(Y3C) real(Y4C) imag(Y4C) real(Y5C) imag(Y5C)];
end



